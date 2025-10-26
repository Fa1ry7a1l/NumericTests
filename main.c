// band_bench_signed_sep.c
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <limits.h>

// ---------------- Параметры ----------------
static const int    n               = 100000;  // размер матрицы
static const int    bw_lower        = 4;       // диагоналей ниже главной
static const int    bw_upper        = 4;       // диагоналей выше главной
static const int    iterationsCount = 200;     // итераций x := A*x + b
static const int    epochs          = 5;       // усреднение
static const int    S               = 32767;   // масштаб для int16_t (|val|<1)

// ---------------- Метаданные ленты ----------------
typedef struct {
    int n, k;
    int *offsets;
    int *start_i;
    int *len;
    size_t *ptr;
} BandMeta;

static void meta_free(BandMeta *m){
    if(!m) return;
    free(m->offsets); free(m->start_i); free(m->len); free(m->ptr);
    m->offsets = m->start_i = m->len = NULL; m->ptr = NULL; m->n = m->k = 0;
}

static BandMeta meta_build(int n, int k, const int *offsets){
    BandMeta M = {0};
    M.n = n; M.k = k;
    M.offsets = (int*)malloc(sizeof(int)*k);
    M.start_i = (int*)malloc(sizeof(int)*k);
    M.len     = (int*)malloc(sizeof(int)*k);
    M.ptr     = (size_t*)malloc(sizeof(size_t)*k);
    size_t total=0;
    for(int d=0; d<k; ++d){
        M.offsets[d] = offsets[d];
        int off = offsets[d];
        int abs_off = (off<0)? -off: off;
        int ld = n - abs_off; if(ld<0) ld=0;
        M.len[d] = ld;
        M.start_i[d] = (off<0)? -off: 0;
        M.ptr[d] = total;
        total += (size_t)ld;
    }
    return M;
}

// ---------------- Форматы хранения ----------------
typedef struct { BandMeta meta; double  *data; }  BandD;
typedef struct { BandMeta meta; int16_t *data; }  BandI16;

// ---------------- Создание/освобождение ----------------
static BandD* bandd_create_random_signed(int n, int bw_lower, int bw_upper){
    int k = bw_lower + 1 + bw_upper;
    int *offs = (int*)malloc(sizeof(int)*k);
    int idx=0; for(int o=-bw_lower; o<=bw_upper; ++o) offs[idx++]=o;

    BandD *A = (BandD*)calloc(1,sizeof(BandD));
    A->meta = meta_build(n, k, offs);
    free(offs);

    size_t total=0; for(int d=0; d<k; ++d) total += (size_t)A->meta.len[d];
    A->data = (double*)malloc(sizeof(double)*total);

    // Случайные значения в (-1,1), главную диагональ слегка усилим
    for(int d=0; d<k; ++d){
        int si = A->meta.start_i[d], ld = A->meta.len[d];
        size_t base = A->meta.ptr[d];
        int off = A->meta.offsets[d];
        for(int t=0; t<ld; ++t){
            double r = ((double)rand()/(double)RAND_MAX)*2.0 - 1.0; // (-1,1)
            if(off==0) r = 0.5*r + (r>=0? 0.5: -0.5);
            if(r >= 1.0) r = 0.999999;
            if(r <= -1.0) r = -0.999999;
            A->data[base + (size_t)t] = r;
        }
    }
    return A;
}
static void bandd_free(BandD *A){ if(!A) return; meta_free(&A->meta); free(A->data); free(A); }

static BandI16* bandi16_from_bandd(const BandD *Ad, int S){
    BandI16 *Aq = (BandI16*)calloc(1,sizeof(BandI16));
    Aq->meta = meta_build(Ad->meta.n, Ad->meta.k, Ad->meta.offsets);
    size_t total=0; for(int d=0; d<Aq->meta.k; ++d) total += (size_t)Aq->meta.len[d];
    Aq->data = (int16_t*)malloc(sizeof(int16_t)*total);

    for(int d=0; d<Aq->meta.k; ++d){
        size_t base = Aq->meta.ptr[d], based = Ad->meta.ptr[d];
        int ld = Aq->meta.len[d];
        for(int t=0; t<ld; ++t){
            long v = lround(Ad->data[based + (size_t)t] * (double)S);
            if(v < -(long)S) v = -(long)S;
            if(v >  (long)S) v =  (long)S;
            Aq->data[base + (size_t)t] = (int16_t)v;
        }
    }
    return Aq;
}
static void bandi16_free(BandI16 *A){ if(!A) return; meta_free(&A->meta); free(A->data); free(A); }

// ---------------- Умножения ----------------
static void band_matvec_double(const BandD *A, const double *x, const double *b, double *y){
    const BandMeta *m=&A->meta; int n=m->n;
    for(int i=0;i<n;++i) y[i]=0.0;
    for(int d=0; d<m->k; ++d){
        int off=m->offsets[d], si=m->start_i[d], ld=m->len[d];
        const double *diag=A->data + m->ptr[d];
        for(int t=0; t<ld; ++t){
            int i=si+t, j=i+off;
            y[i] += diag[t] * x[j];
        }
    }
    if(b) for(int i=0;i<n;++i) y[i]+=b[i];
}

// Чисто целочисленное умножение: yS2 = A_S * x_S + b_S2, где
//   A_S = round(A*S), x_S = round(x*S), b_S2 = round(b*S^2).
// НИЧЕГО не делим внутри — только целочисленные операции и суммирование.
static void band_matvec_i16_pure(const BandI16 *Aq,
                                 const int16_t *xS,
                                 const int32_t *bS2,
                                 int64_t *yS2_out, // результат в масштабе S^2
                                 int add_bias)     // 1 => добавить bS2, 0 => нет
{
    const BandMeta *m=&Aq->meta; int n=m->n;
    // обнулим
    for(int i=0;i<n;++i) yS2_out[i] = 0;

    for(int d=0; d<m->k; ++d){
        int off=m->offsets[d], si=m->start_i[d], ld=m->len[d];
        const int16_t *diagS=Aq->data + m->ptr[d];
        for(int t=0; t<ld; ++t){
            int i=si+t, j=i+off;
            // prod = aS * xS  (|prod| <= S^2), копим в int64_t
            int32_t prod = (int32_t)diagS[t] * (int32_t)xS[j];
            yS2_out[i] += (int64_t)prod;
        }
    }
    if(add_bias && bS2){
        for(int i=0;i<n;++i) yS2_out[i] += (int64_t)bS2[i];
    }
}

// ---------------- Квантизация x, b и деквантизация ----------------
static void quantize_x_b_i16(const double *x, const double *b, int n,
                             int16_t *xS, int32_t *bS2, int S){
    const long long S2 = (long long)S*(long long)S;
    for(int i=0;i<n;++i){
        long vx = lround(x[i] * (double)S);
        if(vx < -(long)S) vx = -(long)S;
        if(vx >  (long)S) vx =  (long)S;
        xS[i] = (int16_t)vx;

        if(bS2){
            long long vb = llround(b[i] * (double)S2);
            if(vb < (long long)INT32_MIN) vb = (long long)INT32_MIN;
            if(vb > (long long)INT32_MAX) vb = (long long)INT32_MAX;
            bS2[i] = (int32_t)vb;
        }
    }
}

static void dequantize_y_update_xS(const int64_t *yS2, double *y,
                                   int16_t *xS_next, int n, int S){
    const double invS2 = 1.0 / ((double)S * (double)S);
    for(int i=0;i<n;++i){
        double val = (double)yS2[i] * invS2;  // деквантизация
        y[i] = val;
        long vx = lround(val * (double)S);   // квантизация для следующей итерации
        if(vx < -(long)S) vx = -(long)S;
        if(vx >  (long)S) vx =  (long)S;
        xS_next[i] = (int16_t)vx;
    }
}

// ---------------- Метрики ----------------
static double rel_l1(const double *ref, const double *test, int n){
    double num=0, den=0; for(int i=0;i<n;++i){ num += fabs(ref[i]-test[i]); den += fabs(ref[i]); }
    return (den>0)? (100.0*num/den) : 0.0;
}
static double rel_l2(const double *ref, const double *test, int n){
    double num=0, den=0; for(int i=0;i<n;++i){ double d=ref[i]-test[i]; num+=d*d; den+=ref[i]*ref[i]; }
    return (den>0)? (100.0*sqrt(num/den)) : 0.0;
}

// ---------------- Главный тест ----------------
int main(void){
    long long seed = (long long)time(NULL);
    printf("seed: %lld\n", seed);
    srand((unsigned)seed);

    // Матрица A (double)
    BandD *A = bandd_create_random_signed(n, bw_lower, bw_upper);

    // Одноразовая квантизация матрицы в int16_t — замеряем отдельно
    clock_t t_mat_quant = clock();
    BandI16 *Aq = bandi16_from_bandd(A, S);
    t_mat_quant = clock() - t_mat_quant;

    // Векторы/буферы
    double  *x    = (double*)malloc(sizeof(double)*n);
    double  *b    = (double*)malloc(sizeof(double)*n);
    double  *y_d  = (double*)malloc(sizeof(double)*n);
    double  *y_q  = (double*)malloc(sizeof(double)*n);
    double  *tmp  = (double*)malloc(sizeof(double)*n);
    int16_t *xS   = (int16_t*)malloc(sizeof(int16_t)*n);
    int32_t *bS2  = (int32_t*)malloc(sizeof(int32_t)*n);
    int16_t *xS2  = (int16_t*)malloc(sizeof(int16_t)*n);
    int64_t *yS2  = (int64_t*)malloc(sizeof(int64_t)*n);

    // Случайные x,b в (-1,1)
    for(int i=0;i<n;++i){
        double rx = ((double)rand()/(double)RAND_MAX)*2.0 - 1.0;
        double rb = ((double)rand()/(double)RAND_MAX)*2.0 - 1.0;
        if(rx >= 1.0) rx = 0.999999; if(rx <= -1.0) rx = -0.999999;
        if(rb >= 1.0) rb = 0.999999; if(rb <= -1.0) rb = -0.999999;
        x[i] = rx; b[i] = rb;
    }

    // -------- baseline (double) --------
    clock_t t_double = 0;
    for(int e=0; e<epochs; ++e){
        for(int i=0;i<n;++i) tmp[i] = x[i];
        clock_t t0 = clock();
        for(int it=0; it<iterationsCount; ++it){
            band_matvec_double(A, tmp, b, y_d);
            for(int i=0;i<n;++i) tmp[i] = y_d[i];
        }
        t_double += clock() - t0;
    }

    // -------- optimized (int16) --------
    // Время на преобразования (кроме матрицы): квантизация x,b и пер-итерационные
    clock_t t_conv = 0;
    // Время чистого целочисленного матричного умножения
    clock_t t_matmul_int = 0;

    for(int e=0; e<epochs; ++e){
        // (a) квантизация начальных x,b — учитываем в t_conv
        clock_t t0 = clock();
        quantize_x_b_i16(x, b, n, xS, bS2, S);
        t_conv += clock() - t0;

        // (b) итерации: раздельно считаем matmul и convert
        for(int it=0; it<iterationsCount; ++it){
            // матрично-векторное умножение ТОЛЬКО в целых
            t0 = clock();
            band_matvec_i16_pure(Aq, xS, bS2, yS2, /*add_bias=*/1);
            t_matmul_int += clock() - t0;

            // деквантизация yS2 -> y_q и квантизация в xS для следующей итерации
            t0 = clock();
            dequantize_y_update_xS(yS2, y_q, xS2, n, S);
            // xS := xS2
            int64_t bytes = (int64_t)n * sizeof(int16_t);
            // простое копирование (без memcpy, чтобы быть честными к clock()):
            for(int i=0;i<n;++i) xS[i] = xS2[i];
            t_conv += clock() - t0;
        }
    }

    // -------- точность (одна итерация) --------
    band_matvec_double(A, x, b, y_d);
    quantize_x_b_i16(x, b, n, xS, bS2, S);
    band_matvec_i16_pure(Aq, xS, bS2, yS2, 1);
    dequantize_y_update_xS(yS2, y_q, xS2, n, S);

    double err_l1 = rel_l1(y_d, y_q, n);
    double err_l2 = rel_l2(y_d, y_q, n);

    // -------- вывод --------
    double tD   = (double)t_double      / CLOCKS_PER_SEC / epochs;
    double tMM  = (double)t_matmul_int  / CLOCKS_PER_SEC / epochs;
    double tCV  = (double)t_conv        / CLOCKS_PER_SEC / epochs;
    double tAQ  = (double)t_mat_quant   / CLOCKS_PER_SEC; // одноразово

    printf("RelError L1: %.6f%%\n", err_l1);
    printf("RelError L2: %.6f%%\n", err_l2);
    printf("Time (double baseline):           %.6f s/epoch\n", tD);
    printf("Time (int16 matmul only):         %.6f s/epoch\n", tMM);
    printf("Time (quant/dequant per-epoch):   %.6f s/epoch\n", tCV);
    printf("Time (matrix quantization once):  %.6f s (one-off)\n", tAQ);
    if(tMM > 0.0) printf("Speedup (double / int16-matmul):  %.2fx\n", tD / tMM);

    // -------- очистка --------
    free(x); free(b); free(y_d); free(y_q); free(tmp);
    free(xS); free(bS2); free(xS2); free(yS2);
    bandd_free(A); bandi16_free(Aq);
    return 0;
}
