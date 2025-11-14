// band_bench_1e6_Ax_band_and_csr.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <limits.h>

// ================== ПАРАМЕТРЫ ТЕСТА ==================
static const int n = 1000000;           // размер матрицы (n x n)
static const int band_half_width = 31; // половина ширины полосы
static const int band_count = 20;       // число полос
static const int band_centers[20] = {0, -100, +100, -200, +200, -300, +300, -400, +400, -500, +500, -600, + 600, -700, +700, -800 , +800,-900,+900,-1000}; // центры полос (offset = j-i)

static const int iterationsCount = 20; // итерации x := A*x
static const int epochs = 3;           // усреднение по эпохам
static const int S = 4096;            // масштаб для int16_t (|val|<1)

#define PARTIAL_CHUNK      128   // размер блока диагоналей для частичных сумм в band-int16
#define CSR_PARTIAL_CHUNK 128   // размер блока nnz для частичных сумм в csr-int16

// ================== МЕТАДАННЫЕ ЛЕНТЫ ==================
typedef struct {
    int n, k;
    int *offsets;
    int *start_i;
    int *len;
    size_t *ptr;
} BandMeta;

static void meta_free(BandMeta *m) {
    if (!m) return;
    free(m->offsets);
    free(m->start_i);
    free(m->len);
    free(m->ptr);
    m->offsets = m->start_i = m->len = NULL;
    m->ptr = NULL;
    m->n = m->k = 0;
}

static int *build_offsets_from_bands(int *out_k) {
    int per_band = 2 * band_half_width + 1;
    int k = band_count * per_band;
    int *offs = (int *) malloc(sizeof(int) * k);
    if (!offs) {
        fprintf(stderr, "build_offsets_from_bands: alloc failed\n");
        exit(1);
    }
    int idx = 0;
    for (int b = 0; b < band_count; ++b) {
        int c = band_centers[b];
        for (int d = -band_half_width; d <= band_half_width; ++d) {
            int off = c + d;
            if (off <= -n + 1) off = -n + 2;
            if (off >=  n - 1) off =  n - 2;
            offs[idx++] = off;
        }
    }
    *out_k = k;
    return offs;
}

static BandMeta meta_build(int n, int k, const int *offsets) {
    BandMeta M = (BandMeta){0};
    M.n = n;
    M.k = k;
    M.offsets = (int *) malloc(sizeof(int) * k);
    M.start_i = (int *) malloc(sizeof(int) * k);
    M.len     = (int *) malloc(sizeof(int) * k);
    M.ptr     = (size_t *) malloc(sizeof(size_t) * k);
    if (!M.offsets || !M.start_i || !M.len || !M.ptr) {
        fprintf(stderr, "meta_build: alloc failed\n");
        exit(1);
    }
    size_t total = 0;
    for (int d = 0; d < k; ++d) {
        M.offsets[d] = offsets[d];
        int off = offsets[d];
        int abs_off = (off < 0) ? -off : off;
        int ld = n - abs_off;
        if (ld < 0) ld = 0;
        M.len[d] = ld;
        M.start_i[d] = (off < 0) ? -off : 0;
        M.ptr[d] = total;
        total += (size_t) ld;
    }
    return M;
}

// ================== ХРАНЕНИЕ ДВУХ ВАРИАНТОВ ЛЕНТЫ ==================
typedef struct {
    BandMeta meta;
    double *data;   // baseline в double
} BandD;

typedef struct {
    BandMeta meta;
    int16_t *data;  // квантованная матрица
} BandI16;

static BandD *bandd_create_random_signed(int n) {
    int k;
    int *offs = build_offsets_from_bands(&k);

    BandD *A = (BandD *) calloc(1, sizeof(BandD));
    if (!A) {
        fprintf(stderr, "bandd_create_random_signed: alloc BandD failed\n");
        exit(1);
    }

    A->meta = meta_build(n, k, offs);
    free(offs);

    size_t total = 0;
    for (int d = 0; d < k; ++d) total += (size_t) A->meta.len[d];
    A->data = (double *) malloc(sizeof(double) * total);
    if (!A->data) {
        fprintf(stderr, "Alloc failed for BandD (%.2f GB)\n",
                total * 8.0 / (1024.0 * 1024.0 * 1024.0));
        exit(1);
    }

    // Случайные (-1,1), усилим центральные диагонали полос
    for (int d = 0; d < k; ++d) {
        int off  = A->meta.offsets[d];
        int ld   = A->meta.len[d];
        size_t base = A->meta.ptr[d];
        int is_center = 0;
        for (int c = 0; c < band_count; ++c) {
            if (off == band_centers[c]) {
                is_center = 1;
                break;
            }
        }
        for (int t = 0; t < ld; ++t) {
            double r = ((double) rand() / (double) RAND_MAX) * 2.0 - 1.0; // (-1,1)
            if (is_center) r = 0.5 * r + (r >= 0 ? 0.5 : -0.5);
            if (r >= 1.0) r = 0.999999;
            if (r <= -1.0) r = -0.999999;
            A->data[base + (size_t) t] = r;
        }
    }
    return A;
}

static BandI16 *bandi16_from_bandd(const BandD *Ad, int S) {
    BandI16 *Aq = (BandI16 *) calloc(1, sizeof(BandI16));
    if (!Aq) {
        fprintf(stderr, "bandi16_from_bandd: alloc BandI16 failed\n");
        exit(1);
    }

    Aq->meta = meta_build(Ad->meta.n, Ad->meta.k, Ad->meta.offsets);

    size_t total = 0;
    for (int d = 0; d < Aq->meta.k; ++d) total += (size_t) Aq->meta.len[d];
    Aq->data = (int16_t *) malloc(sizeof(int16_t) * total);
    if (!Aq->data) {
        fprintf(stderr, "Alloc failed for BandI16 (%.2f GB)\n",
                total * 2.0 / (1024.0 * 1024.0 * 1024.0));
        exit(1);
    }

    for (int d = 0; d < Aq->meta.k; ++d) {
        size_t base  = Aq->meta.ptr[d];
        size_t based = Ad->meta.ptr[d];
        int ld = Aq->meta.len[d];
        for (int t = 0; t < ld; ++t) {
            double scaled = Ad->data[based + (size_t)t] * (double)S;
            Aq->data[base + (size_t)t] = (int16_t)scaled; // предполагается, что не переполняется
        }
    }
    return Aq;
}

static void bandd_free(BandD *A) {
    if (!A) return;
    meta_free(&A->meta);
    free(A->data);
    free(A);
}

static void bandi16_free(BandI16 *A) {
    if (!A) return;
    meta_free(&A->meta);
    free(A->data);
    free(A);
}

// ================== CSR ДЛЯ int16 ==================
typedef struct {
    int n;
    size_t nnz;
    int *row_ptr;    // размер n+1
    int *col_idx;    // размер nnz
    int16_t *vals;   // размер nnz (int16 матрица)
} CSR_I16;

static void csr_free(CSR_I16 *A) {
    if (!A) return;
    free(A->row_ptr);
    free(A->col_idx);
    free(A->vals);
    free(A);
}

// Строим CSR из ленточного представления int16, пропуская нули.
static CSR_I16 *csr_from_bandI16(const BandI16 *Aq, clock_t *t_build_out) {
    const BandMeta *m = &Aq->meta;
    int n = m->n;
    int k = m->k;

    clock_t t0 = clock();

    int *row_counts = (int *) calloc((size_t)n, sizeof(int));
    if (!row_counts) {
        fprintf(stderr, "csr_from_bandI16: alloc row_counts failed\n");
        exit(1);
    }

    long long nnz_ll = 0;
    for (int d = 0; d < k; ++d) {
        int off  = m->offsets[d];
        int si   = m->start_i[d];
        int ld   = m->len[d];
        size_t base = m->ptr[d];
        const int16_t *diagS = Aq->data + base;
        for (int t = 0; t < ld; ++t) {
            int16_t val = diagS[t];
            if (val == 0) continue;
            int i = si + t;
            row_counts[i]++;
            nnz_ll++;
        }
    }

    if (nnz_ll < 0) {
        fprintf(stderr, "csr_from_bandI16: nnz overflow\n");
        exit(1);
    }
    size_t nnz = (size_t)nnz_ll;

    CSR_I16 *C = (CSR_I16 *) calloc(1, sizeof(CSR_I16));
    if (!C) {
        fprintf(stderr, "csr_from_bandI16: alloc CSR_I16 failed\n");
        exit(1);
    }
    C->n   = n;
    C->nnz = nnz;
    C->row_ptr = (int *) malloc(sizeof(int) * (size_t)(n+1));
    C->col_idx = (int *) malloc(sizeof(int) * nnz);
    C->vals    = (int16_t *) malloc(sizeof(int16_t) * nnz);
    if (!C->row_ptr || !C->col_idx || !C->vals) {
        fprintf(stderr, "csr_from_bandI16: alloc CSR arrays failed\n");
        exit(1);
    }

    C->row_ptr[0] = 0;
    for (int i = 0; i < n; ++i) {
        C->row_ptr[i+1] = C->row_ptr[i] + row_counts[i];
    }

    int *row_pos = (int *) malloc(sizeof(int) * (size_t)n);
    if (!row_pos) {
        fprintf(stderr, "csr_from_bandI16: alloc row_pos failed\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) {
        row_pos[i] = C->row_ptr[i];
    }

    for (int d = 0; d < k; ++d) {
        int off  = m->offsets[d];
        int si   = m->start_i[d];
        int ld   = m->len[d];
        size_t base = m->ptr[d];
        const int16_t *diagS = Aq->data + base;
        for (int t = 0; t < ld; ++t) {
            int16_t val = diagS[t];
            if (val == 0) continue;
            int i = si + t;
            int j = i + off;
            int pos = row_pos[i]++;
            C->col_idx[pos] = j;
            C->vals[pos]    = val;
        }
    }

    free(row_counts);
    free(row_pos);

    clock_t t1 = clock();
    if (t_build_out) *t_build_out = (t1 - t0);

    return C;
}

// ================== УМНОЖЕНИЕ ==================
// 1) baseline: y = A(double) * x
static void band_matvec_double(const BandD *A, const double *x, double *y) {
    const BandMeta *m = &A->meta;
    int n = m->n;
    for (int i = 0; i < n; ++i) y[i] = 0.0;
    for (int d = 0; d < m->k; ++d) {
        int off = m->offsets[d], si = m->start_i[d], ld = m->len[d];
        const double *diag = A->data + m->ptr[d];
        for (int t = 0; t < ld; ++t) {
            int i = si + t, j = i + off;
            y[i] += diag[t] * x[j];
        }
    }
}

// 2) optimized band (int16): yS2 (int32) = A_S * x_S, с частичными суммами
static void band_matvec_i16_pure(const BandI16 *Aq,
                                 const int16_t *xS,
                                 int32_t *yS2_out) {
    const BandMeta *m = &Aq->meta;
    int n = m->n;
    int k = m->k;

    const int max_slots = (k + PARTIAL_CHUNK - 1) / PARTIAL_CHUNK;

    for (int i = 0; i < n; ++i) {
        int32_t t[max_slots];
        int t_count = 0;

        for (int base = 0; base < k; base += PARTIAL_CHUNK) {
            int upto = base + PARTIAL_CHUNK;
            if (upto > k) upto = k;

            int32_t local = 0;
            for (int d = base; d < upto; ++d) {
                int si  = m->start_i[d];
                int ld  = m->len[d];
                if (i < si || i >= si + ld) continue;

                int off = m->offsets[d];
                int ti  = i - si;
                int j   = i + off;

                const int16_t *diagS = Aq->data + m->ptr[d];
                int32_t a = (int32_t)diagS[ti];
                int32_t b = (int32_t)xS[j];
                local += a * b;   // предполагается, что не переполняется
            }
            t[t_count++] = local; // временная частичная сумма в int
        }

        int32_t acc = 0;
        for (int s = 0; s < t_count; ++s) {
            acc += t[s];          // финальная сумма в int
        }
        yS2_out[i] = acc;
    }
}

// 3) optimized CSR (int16): yS2 (int32) = A_S * x_S по CSR
static void csr_matvec_i16(const CSR_I16 *A,
                           const int16_t *xS,
                           int32_t *yS2_out) {
    int n = A->n;
    const int *row_ptr = A->row_ptr;
    const int *col_idx = A->col_idx;
    const int16_t *vals = A->vals;

    for (int i = 0; i < n; ++i) {
        int row_start = row_ptr[i];
        int row_end   = row_ptr[i+1];
        int nnz_row   = row_end - row_start;

        int max_slots = (nnz_row + CSR_PARTIAL_CHUNK - 1) / CSR_PARTIAL_CHUNK;
        if (max_slots <= 0) {
            yS2_out[i] = 0;
            continue;
        }

        int32_t t[max_slots];
        int t_count = 0;

        int p = row_start;
        while (p < row_end) {
            int32_t local = 0;
            int limit = p + CSR_PARTIAL_CHUNK;
            if (limit > row_end) limit = row_end;

            for (; p < limit; ++p) {
                int j = col_idx[p];
                int32_t a = (int32_t)vals[p];
                int32_t b = (int32_t)xS[j];
                local += a * b;    // предполагается, что не переполняется
            }

            t[t_count++] = local;
        }

        int32_t acc = 0;
        for (int s = 0; s < t_count; ++s) {
            acc += t[s];
        }
        yS2_out[i] = acc;
    }
}

// ================== КОНВЕРСИИ (время считаем отдельно) ==================
static void quantize_x_i16(const double *x, int n,
                           int16_t *xS, int S) {
    (void)S;
    for (int i = 0; i < n; ++i) {
        double scaled = x[i] * (double)S;
        xS[i] = (int16_t)scaled; // предполагается, что не переполняется
    }
}

static void dequantize_y_update_xS(const int32_t *yS2, double *y,
                                   int16_t *xS_next, int n, int S) {
    const double invS2 = 1.0 / ((double) S * (double) S);
    for (int i = 0; i < n; ++i) {
        double val = (double)yS2[i] * invS2;
        y[i] = val;
        double scaled = val * (double)S;
        xS_next[i] = (int16_t)scaled; // предполагается, что не переполняется
    }
}

// ================== МЕТРИКИ ==================
static double rel_l1(const double *ref, const double *test, int n) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
        num += fabs(ref[i] - test[i]);
        den += fabs(ref[i]);
    }
    return (den > 0.0) ? (100.0 * num / den) : 0.0;
}

static double rel_l2(const double *ref, const double *test, int n) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = ref[i] - test[i];
        num += d * d;
        den += ref[i] * ref[i];
    }
    return (den > 0.0) ? (100.0 * sqrt(num / den)) : 0.0;
}

// ================== MAIN ==================
int main(void) {
    long long seed = (long long) time(NULL);
    printf("seed: %lld\n", seed);
    srand((unsigned) seed);

    // 1) Матрица в double
    BandD *A = bandd_create_random_signed(n);

    // 2) Квантование матрицы в int16 (лента)
    clock_t t_mat_quant = clock();
    BandI16 *Aq = bandi16_from_bandd(A, S);
    t_mat_quant = clock() - t_mat_quant;

    // 3) Построение CSR из int16-ленты
    clock_t t_csr_build = 0;
    CSR_I16 *Ac = csr_from_bandI16(Aq, &t_csr_build);

    // Векторы/буферы
    double *x           = (double *) malloc(sizeof(double) * n);
    double *y_d         = (double *) malloc(sizeof(double) * n); // baseline
    double *y_q_band    = (double *) malloc(sizeof(double) * n); // band-int16 деквант
    double *y_q_csr     = (double *) malloc(sizeof(double) * n); // csr-int16 деквант
    double *tmp         = (double *) malloc(sizeof(double) * n);

    int16_t *xS_band    = (int16_t *) malloc(sizeof(int16_t) * n);
    int16_t *xS2_band   = (int16_t *) malloc(sizeof(int16_t) * n);
    int32_t *yS2_band   = (int32_t *) malloc(sizeof(int32_t) * n);

    int16_t *xS_csr     = (int16_t *) malloc(sizeof(int16_t) * n);
    int16_t *xS2_csr    = (int16_t *) malloc(sizeof(int16_t) * n);
    int32_t *yS2_csr    = (int32_t *) malloc(sizeof(int32_t) * n);

    if (!x || !y_d || !y_q_band || !y_q_csr || !tmp ||
        !xS_band || !xS2_band || !yS2_band ||
        !xS_csr  || !xS2_csr  || !yS2_csr) {
        fprintf(stderr, "Allocation failed\n");
        return 1;
    }

    // Случайный x в (-1,1)
    for (int i = 0; i < n; ++i) {
        double rx = ((double) rand() / (double) RAND_MAX) * 2.0 - 1.0;
        if (rx >= 1.0) rx = 0.999999;
        if (rx <= -1.0) rx = -0.999999;
        x[i] = rx;
    }

    // ----- 1) baseline (double): A*x -----
    clock_t t_double = 0;
    for (int e = 0; e < epochs; ++e) {
        for (int i = 0; i < n; ++i) tmp[i] = x[i];
        clock_t t0 = clock();
        for (int it = 0; it < iterationsCount; ++it) {
            band_matvec_double(A, tmp, y_d);
            for (int i = 0; i < n; ++i) tmp[i] = y_d[i];
        }
        t_double += clock() - t0;
    }

    // ----- 2) band-int16: A*x -----
    clock_t t_conv_band   = 0;
    clock_t t_matmul_band = 0;

    for (int e = 0; e < epochs; ++e) {
        clock_t t0 = clock();
        quantize_x_i16(x, n, xS_band, S);
        t_conv_band += clock() - t0;

        for (int it = 0; it < iterationsCount; ++it) {
            t0 = clock();
            band_matvec_i16_pure(Aq, xS_band, yS2_band);
            t_matmul_band += clock() - t0;

            t0 = clock();
            dequantize_y_update_xS(yS2_band, y_q_band, xS2_band, n, S);
            for (int i = 0; i < n; ++i) xS_band[i] = xS2_band[i];
            t_conv_band += clock() - t0;
        }
    }

    // ----- 3) CSR-int16: A*x -----
    clock_t t_conv_csr   = 0;
    clock_t t_matmul_csr = 0;

    for (int e = 0; e < epochs; ++e) {
        clock_t t0 = clock();
        quantize_x_i16(x, n, xS_csr, S);
        t_conv_csr += clock() - t0;

        for (int it = 0; it < iterationsCount; ++it) {
            t0 = clock();
            csr_matvec_i16(Ac, xS_csr, yS2_csr);
            t_matmul_csr += clock() - t0;

            t0 = clock();
            dequantize_y_update_xS(yS2_csr, y_q_csr, xS2_csr, n, S);
            for (int i = 0; i < n; ++i) xS_csr[i] = xS2_csr[i];
            t_conv_csr += clock() - t0;
        }
    }

    // ----- Оценка точности на одной итерации (y = A*x) -----
    band_matvec_double(A, x, y_d);

    quantize_x_i16(x, n, xS_band, S);
    band_matvec_i16_pure(Aq, xS_band, yS2_band);
    dequantize_y_update_xS(yS2_band, y_q_band, xS2_band, n, S);

    quantize_x_i16(x, n, xS_csr, S);
    csr_matvec_i16(Ac, xS_csr, yS2_csr);
    dequantize_y_update_xS(yS2_csr, y_q_csr, xS2_csr, n, S);

    double err_l1_band = rel_l1(y_d, y_q_band, n);
    double err_l2_band = rel_l2(y_d, y_q_band, n);
    double err_l1_csr  = rel_l1(y_d, y_q_csr,  n);
    double err_l2_csr  = rel_l2(y_d, y_q_csr,  n);

    // ----- вывод таймингов и ошибок -----
    double tD    = (double) t_double      / CLOCKS_PER_SEC / epochs;
    double tMM_b = (double) t_matmul_band / CLOCKS_PER_SEC / epochs;
    double tCV_b = (double) t_conv_band   / CLOCKS_PER_SEC / epochs;
    double tMM_c = (double) t_matmul_csr  / CLOCKS_PER_SEC / epochs;
    double tCV_c = (double) t_conv_csr    / CLOCKS_PER_SEC / epochs;
    double tAQ   = (double) t_mat_quant   / CLOCKS_PER_SEC;
    double tCSR  = (double) t_csr_build   / CLOCKS_PER_SEC;

    printf("n = %d, bands = %d x ~%d diagonals\n", n, band_count, 2 * band_half_width + 1);
    printf("S = %d\n",S);
    printf("PARTIAL_CHUNK (band-int16) = %d\n", PARTIAL_CHUNK);
    printf("CSR_PARTIAL_CHUNK          = %d\n", CSR_PARTIAL_CHUNK);

    printf("RelError (band-int16) L1: %.6f%%\n", err_l1_band);
    printf("RelError (band-int16) L2: %.6f%%\n", err_l2_band);
    printf("RelError (csr-int16 ) L1: %.6f%%\n", err_l1_csr);
    printf("RelError (csr-int16 ) L2: %.6f%%\n", err_l2_csr);

    printf("Time (double baseline):                 %.6f s/epoch\n", tD);
    printf("Time (band-int16 matmul only):          %.6f s/epoch\n", tMM_b);
    printf("Time (band-int16 quant/dequant epoch):  %.6f s/epoch\n", tCV_b);
    printf("Time (csr-int16 matmul only):           %.6f s/epoch\n", tMM_c);
    printf("Time (csr-int16 quant/dequant epoch):   %.6f s/epoch\n", tCV_c);
    printf("Time (matrix quantization once):        %.6f s\n",      tAQ);
    printf("Time (CSR build from band-int16 once):  %.6f s\n",      tCSR);

    if (tMM_b > 0.0) printf("Speedup (double / band-int16-matmul):  %.2fx\n", tD / tMM_b);
    if (tMM_c > 0.0) printf("Speedup (double / csr-int16-matmul):   %.2fx\n", tD / tMM_c);

    // ----- Суммы результатов (для проверки, что вектора "живые") -----
    double sum_y_d      = 0.0;
    double sum_y_band   = 0.0;
    double sum_y_csr    = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_y_d    += y_d[i];
        sum_y_band += y_q_band[i];
        sum_y_csr  += y_q_csr[i];
    }
    printf("Sum(y_double baseline):     %.10e\n", sum_y_d);
    printf("Sum(y_band int16 dequant):  %.10e\n", sum_y_band);
    printf("Sum(y_csr  int16 dequant):  %.10e\n", sum_y_csr);

    // ----- очистка -----
    free(x);
    free(y_d);
    free(y_q_band);
    free(y_q_csr);
    free(tmp);

    free(xS_band);
    free(xS2_band);
    free(yS2_band);

    free(xS_csr);
    free(xS2_csr);
    free(yS2_csr);

    bandd_free(A);
    bandi16_free(Aq);
    csr_free(Ac);

    return 0;
}
