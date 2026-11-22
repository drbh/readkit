/*
 * Minimal C ABI over the NCBI VDB C API, compiled only with the `ngs`
 * feature. The VDB symbols are exported by the SRA Toolkit's
 * `libncbi-ngs` library; the
 * toolkit does not install VDB headers, so the few prototypes used here
 * are declared below, copied from the ncbi-vdb headers under `interfaces/vdb/`.
 *
 * Why VDB and not the NGS C++ `ReadIterator`: NGS
 * fragments are biological reads only and expose no read type, so
 * technical reads (single-cell barcodes) and orientation bits are
 * unreachable through NGS.
 *
 * One handle = one SEQUENCE table, one read cursor. Handles are not
 * thread-safe; use one per thread. All handles share one VDB manager, made
 * once and never released (see `shared_manager`). Pointers returned by
 * `rk_read` point into VDB's cursor cache and stay valid until the next
 * `rk_read` or `rk_close` on the same handle.
 */

#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef uint32_t rc_t;
typedef struct VDBManager VDBManager;
typedef struct VDatabase VDatabase;
typedef struct VTable VTable;
typedef struct VCursor VCursor;
typedef struct VSchema VSchema;
typedef struct KDirectory KDirectory;

rc_t VDBManagerMakeRead(const VDBManager **mgr, const KDirectory *wd);
rc_t VDBManagerRelease(const VDBManager *self);
rc_t VDBManagerGetQualityString(const VDBManager *self, const char **quality);
rc_t VDBManagerOpenDBRead(const VDBManager *self, const VDatabase **db,
                          const VSchema *schema, const char *path, ...);
rc_t VDBManagerOpenTableRead(const VDBManager *self, const VTable **tbl,
                             const VSchema *schema, const char *path, ...);
rc_t VDatabaseOpenTableRead(const VDatabase *self, const VTable **tbl,
                            const char *name, ...);
rc_t VDatabaseRelease(const VDatabase *self);
rc_t VDatabaseGetQualityCapability(const VDatabase *self, bool *full, bool *synth);
rc_t VTableGetQualityCapability(const VTable *self, bool *full, bool *synth);
rc_t VTableRelease(const VTable *self);
rc_t VTableCreateCursorRead(const VTable *self, const VCursor **curs);
rc_t VCursorAddColumn(const VCursor *self, uint32_t *idx, const char *name, ...);
rc_t VCursorOpen(const VCursor *self);
rc_t VCursorIdRange(const VCursor *self, uint32_t idx, int64_t *first, uint64_t *count);
rc_t VCursorCellDataDirect(const VCursor *self, int64_t row_id, uint32_t col_idx,
                           uint32_t *elem_bits, const void **base, uint32_t *boff,
                           uint32_t *row_len);
rc_t VCursorRelease(const VCursor *self);
rc_t string_printf(char *dst, size_t bsize, size_t *num_writ, const char *fmt, ...);
rc_t VTableTypespec(const VTable *self, char *ts_buff, size_t ts_size);

/* KDB: raw physical blobs, for the fast path in src/fast.rs. */
typedef struct KTable KTable;
typedef struct KColumn KColumn;
typedef struct KColumnBlob KColumnBlob;
rc_t VTableOpenKTableRead(const VTable *self, const KTable **ktbl);
rc_t KTableRelease(const KTable *self);
rc_t KTableOpenColumnRead(const KTable *self, const KColumn **col, const char *name, ...);
rc_t KColumnRelease(const KColumn *self);
rc_t KColumnOpenBlobRead(const KColumn *self, const KColumnBlob **blob, int64_t id);
rc_t KColumnBlobIdRange(const KColumnBlob *self, int64_t *first, uint32_t *count);
rc_t KColumnBlobRead(const KColumnBlob *self, size_t offset, void *buffer, size_t bsize,
                     size_t *num_read, size_t *remaining);
rc_t KColumnBlobRelease(const KColumnBlob *self);

/* Mirrors `ffi::Spot` in src/ngs.rs. */
typedef struct rk_spot {
    const char *bases;
    const char *quals;
    const int32_t *read_start;
    const uint32_t *read_len;
    const uint8_t *read_type;
    const char *name;
    uint32_t bases_len;
    uint32_t quals_len;
    uint32_t nreads;
    uint32_t nstarts;
    uint32_t nlens;
    uint32_t name_len;
} rk_spot;

/* Mirrors `ffi::Info` in src/ngs.rs. */
typedef struct rk_info {
    int64_t first_row;
    uint64_t row_count;
    /* 0 = unknown, 1 = original (full) qualities, 2 = synthetic only (SRA Lite) */
    int32_t quality;
    /* 1 = the run is a database (cSRA or similar), 0 = a bare table */
    int32_t is_database;
    int32_t has_names;
} rk_info;

/* Mirrors `ffi::Status` in src/ngs.rs. */
typedef struct rk_status {
    rc_t rc;
    /* which step failed, a static string */
    const char *step;
    char message[512];
} rk_status;

enum { C_READ, C_QUAL, C_START, C_LEN, C_TYPE, C_NAME, C_COUNT };

typedef struct rk_handle {
    const VDBManager *mgr;
    const VDatabase *db;
    const VTable *tbl;
    const VCursor *curs;
    uint32_t col[C_COUNT];
    /* 0: READ and QUALITY are not in the cursor (the fast path reads them
     * from KDB blobs). */
    int want_seq;
    rk_info info;
} rk_handle;

static int fail(rk_status *st, rc_t rc, const char *step) {
    if (st) {
        size_t n = 0;
        st->rc = rc;
        st->step = step;
        st->message[0] = '\0';
        if (rc != 0 && string_printf(st->message, sizeof st->message, &n, "%R", rc) != 0) {
            st->message[0] = '\0';
        }
        st->message[sizeof st->message - 1] = '\0';
    }
    return -1;
}

/* One VDB manager for the process, made once and never released.
 *
 * ncbi-vdb 3.2 keeps a process-global quality-preference string that every
 * VDBManagerWhack frees (VDBManagerWhackStatic) and the first database open
 * lazily allocates, neither under a lock. Separate managers made and
 * released on different threads therefore double-free it (seen as
 * "pointer being freed was not allocated" in VDBManagerWhack when readers
 * on several threads close together). A single manager that is never
 * whacked avoids both races; the string is loaded here, under the once. */
static pthread_once_t g_mgr_once = PTHREAD_ONCE_INIT;
static const VDBManager *g_mgr = NULL;
static rc_t g_mgr_rc = 0;

static void make_shared_manager(void) {
    g_mgr_rc = VDBManagerMakeRead(&g_mgr, NULL);
    if (g_mgr_rc == 0) {
        const char *quality = NULL;
        (void)VDBManagerGetQualityString(g_mgr, &quality);
    }
}

static rc_t shared_manager(const VDBManager **mgr) {
    pthread_once(&g_mgr_once, make_shared_manager);
    *mgr = g_mgr;
    return g_mgr_rc;
}

void rk_close(void *handle) {
    rk_handle *h = (rk_handle *)handle;
    if (!h) return;
    if (h->curs) VCursorRelease(h->curs);
    if (h->tbl) VTableRelease(h->tbl);
    if (h->db) VDatabaseRelease(h->db);
    /* h->mgr is the shared manager: not released. */
    free(h);
}

void *rk_open(const char *path, int want_names, int want_seq, rk_info *info, rk_status *st) {
    rk_handle *h = (rk_handle *)calloc(1, sizeof *h);
    rc_t rc;
    bool full = false, synth = false;
    if (!h) {
        fail(st, 0, "allocating handle");
        return NULL;
    }
    h->want_seq = want_seq;
    rc = shared_manager(&h->mgr);
    if (rc) goto err_mgr;
    /* A run is either a database with a SEQUENCE table (cSRA, most newer
     * runs) or a bare SEQUENCE table (older SRA). */
    if (VDBManagerOpenDBRead(h->mgr, &h->db, NULL, "%s", path) == 0) {
        rc = VDatabaseOpenTableRead(h->db, &h->tbl, "%s", "SEQUENCE");
        if (rc) goto err_tbl;
        h->info.is_database = 1;
        rc = VDatabaseGetQualityCapability(h->db, &full, &synth);
    } else {
        h->db = NULL;
        rc = VDBManagerOpenTableRead(h->mgr, &h->tbl, NULL, "%s", path);
        if (rc) goto err_tbl;
        rc = VTableGetQualityCapability(h->tbl, &full, &synth);
    }
    h->info.quality = rc ? 0 : (full ? 1 : (synth ? 2 : 0));

    rc = VTableCreateCursorRead(h->tbl, &h->curs);
    if (rc) goto err_cursor;
    {
        static const char *const cols[C_NAME] = {
            "(INSDC:dna:text)READ",
            "(INSDC:quality:text:phred_33)QUALITY",
            "(INSDC:coord:zero)READ_START",
            "(INSDC:coord:len)READ_LEN",
            "(INSDC:SRA:xread_type)READ_TYPE",
        };
        int i;
        for (i = want_seq ? 0 : C_START; i < C_NAME; i++) {
            rc = VCursorAddColumn(h->curs, &h->col[i], "%s", cols[i]);
            if (rc) {
                fail(st, rc, cols[i]);
                rk_close(h);
                return NULL;
            }
        }
    }
    if (want_names) {
        /* NAME is optional: many runs store no original names. */
        h->info.has_names = VCursorAddColumn(h->curs, &h->col[C_NAME], "%s", "(ascii)NAME") == 0;
    }
    rc = VCursorOpen(h->curs);
    if (rc) goto err_cursor;
    rc = VCursorIdRange(h->curs, 0, &h->info.first_row, &h->info.row_count);
    if (rc) {
        fail(st, rc, "reading the row range");
        rk_close(h);
        return NULL;
    }
    *info = h->info;
    return h;

err_mgr:
    fail(st, rc, "creating the VDB manager");
    rk_close(h);
    return NULL;
err_tbl:
    fail(st, rc, "opening the SEQUENCE table");
    rk_close(h);
    return NULL;
err_cursor:
    fail(st, rc, "opening a cursor");
    rk_close(h);
    return NULL;
}

static int cell(const rk_handle *h, int64_t row, int col, uint32_t bits, const void **base,
                uint32_t *len, rk_status *st, const char *what) {
    uint32_t elem_bits = 0, boff = 0;
    rc_t rc = VCursorCellDataDirect(h->curs, row, h->col[col], &elem_bits, base, &boff, len);
    if (rc) return fail(st, rc, what);
    if (elem_bits != bits || boff != 0) return fail(st, 0, what);
    return 0;
}

/* 0 on success, -1 on error (details in *st). */
int rk_read(void *handle, int64_t row, rk_spot *out, rk_status *st) {
    rk_handle *h = (rk_handle *)handle;
    const void *p;
    uint32_t n;
    memset(out, 0, sizeof *out);
    if (h->want_seq) {
        if (cell(h, row, C_READ, 8, &p, &n, st, "READ")) return -1;
        out->bases = (const char *)p;
        out->bases_len = n;
        if (cell(h, row, C_QUAL, 8, &p, &n, st, "QUALITY")) return -1;
        out->quals = (const char *)p;
        out->quals_len = n;
    }
    if (cell(h, row, C_START, 32, &p, &n, st, "READ_START")) return -1;
    out->read_start = (const int32_t *)p;
    out->nstarts = n;
    if (cell(h, row, C_LEN, 32, &p, &n, st, "READ_LEN")) return -1;
    out->read_len = (const uint32_t *)p;
    out->nlens = n;
    if (cell(h, row, C_TYPE, 8, &p, &n, st, "READ_TYPE")) return -1;
    out->read_type = (const uint8_t *)p;
    out->nreads = n;
    if (h->info.has_names) {
        if (cell(h, row, C_NAME, 8, &p, &n, st, "NAME")) return -1;
        out->name = (const char *)p;
        out->name_len = n;
    }
    return 0;
}

/* The SEQUENCE table's type, e.g. "NCBI:align:tbl:seq#1.1". 0 on success. */
int rk_typespec(void *handle, char *buf, size_t size, rk_status *st) {
    rk_handle *h = (rk_handle *)handle;
    rc_t rc = VTableTypespec(h->tbl, buf, size);
    return rc ? fail(st, rc, "reading the table type") : 0;
}

typedef struct rk_kcol {
    const KTable *ktbl;
    const KColumn *kcol;
} rk_kcol;

void rk_kcol_close(void *kc) {
    rk_kcol *k = (rk_kcol *)kc;
    if (!k) return;
    if (k->kcol) KColumnRelease(k->kcol);
    if (k->ktbl) KTableRelease(k->ktbl);
    free(k);
}

/* The physical column `name` of the SEQUENCE table, or NULL (rc in *st) if
 * the table has none. */
void *rk_kcol_open(void *handle, const char *name, rk_status *st) {
    rk_handle *h = (rk_handle *)handle;
    rk_kcol *k = (rk_kcol *)calloc(1, sizeof *k);
    rc_t rc;
    if (!k) {
        fail(st, 0, "allocating a column handle");
        return NULL;
    }
    rc = VTableOpenKTableRead(h->tbl, &k->ktbl);
    if (rc == 0) rc = KTableOpenColumnRead(k->ktbl, &k->kcol, "%s", name);
    if (rc) {
        fail(st, rc, "opening a physical column");
        rk_kcol_close(k);
        return NULL;
    }
    return k;
}

/* The raw bytes of the blob holding `row`: its id range in first and count and
 * its size in size. Returns 0 when the blob was copied into `buf`, 1 when it
 * is larger than `cap` (nothing copied; call again with *size bytes), -1 on
 * error. */
int rk_kcol_blob(void *kc, int64_t row, int64_t *first, uint32_t *count, uint8_t *buf,
                 size_t cap, size_t *size, rk_status *st) {
    rk_kcol *k = (rk_kcol *)kc;
    const KColumnBlob *b = NULL;
    size_t n = 0, rem = 0, at = 0;
    rc_t rc = KColumnOpenBlobRead(k->kcol, &b, row);
    if (rc) return fail(st, rc, "opening a blob");
    rc = KColumnBlobIdRange(b, first, count);
    if (rc == 0) rc = KColumnBlobRead(b, 0, NULL, 0, &n, &rem);
    if (rc) {
        KColumnBlobRelease(b);
        return fail(st, rc, "reading a blob");
    }
    *size = n + rem;
    if (*size > cap) {
        KColumnBlobRelease(b);
        return 1;
    }
    while (at < *size) {
        rc = KColumnBlobRead(b, at, buf + at, *size - at, &n, &rem);
        if (rc || n == 0) {
            KColumnBlobRelease(b);
            return fail(st, rc, "reading a blob");
        }
        at += n;
    }
    KColumnBlobRelease(b);
    return 0;
}
