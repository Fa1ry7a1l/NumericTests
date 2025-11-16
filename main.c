// perim_only_bench.c
// Тесты только для perim_matvec_double_ref и perim_matvec_i16_ref.

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// ================== ПАРАМЕТРЫ ТЕСТА ==================
static const int n = 1000000;             // размер матрицы (n x n)

// задаём "квазиполосы" через набор диагоналей (offset = j-i)
static const int band_half_width = 31;
static const int band_count = 20;
static const int band_centers[20] = {
    0, -100, +100, -200, +200, -300, +300, -400, +400, -500,
    +500, -600, +600, -700, +700, -800, +800, -900, +900, -1000
};

static const int iterationsCount = 20;    // итерации x := A*x
static const int epochs = 3;              // усреднение
static const int S = 4096;                // масштаб для int16 (|val|<1, * S)

// ================== ФОРМАТ "ПЕРИМЕТР" ==================
// data: все ненулевые элементы подряд (диагональ за диагональю)
// meta: [delta_up, width] * strips
// обход периметра: снизу слева вверх (0..n-1), затем вправо по верхнему краю (n..2n-2)

typedef struct {
    int n;
    int strips;   // число "полос" (диагоналей)
    int *meta;    // длина = 2*strips: [delta_up, width]
    double *data; // все значения double, подряд
} PerimD;

typedef struct {
    int n;
    int strips;
    int *meta;     // тот же meta
    int16_t *data; // значения int16
} PerimI16;

static void perimd_free(PerimD *P){
    if(!P) return;
    free(P->meta);
    free(P->data);
    free(P);
}
static void perimi16_free(PerimI16 *P){
    if(!P) return;
    free(P->meta);
    free(P->data);
    free(P);
}

// ================== ВСПОМОГАТЕЛЬНЫЕ УТИЛИТЫ ==================
static inline int perimeter_contact_pos(int n, int off){
    // путь: левый край i=n-1..0 => позиции 0..n-1
    // затем верхняя строка j=1..n-1 => позиции n..2n-2
    if(off <= 0) {
        // контакт на левом краю: i = -off
        return (n - 1) + off; // в [0..n-1]
    } else {
        // контакт на верхнем краю: j = off
        return n + off;       // в [n..2n-2]
    }
}

static double rel_l1(const double *ref, const double *test, int n){
    double num=0.0, den=0.0;
    for(int i=0;i<n;++i){
        num += fabs(ref[i]-test[i]);
        den += fabs(ref[i]);
    }
    return (den>0.0)?(100.0*num/den):0.0;
}
static double rel_l2(const double *ref, const double *test, int n){
    double num=0.0, den=0.0;
    for(int i=0;i<n;++i){
        double d = ref[i]-test[i];
        num += d*d;
        den += ref[i]*ref[i];
    }
    return (den>0.0)?(100.0*sqrt(num/den)):0.0;
}

static void quantize_x_i16(const double *x, int n, int16_t *xS, int S){
    (void)S;
    for(int i=0;i<n;++i){
        double scaled = x[i]*(double)S;
        xS[i] = (int16_t)scaled; // предполагается, что не переполняется
    }
}

static void dequantize_y_update_xS(const int32_t *yS2, double *y,
                                   int16_t *xS_next, int n, int S){
    double invS2 = 1.0 / ((double)S * (double)S);
    for(int i=0;i<n;++i){
        double val = (double)yS2[i] * invS2;
        y[i] = val;
        double scaled = val*(double)S;
        xS_next[i] = (int16_t)scaled;
    }
}

// ================== ПОСТРОИТЕЛИ PERIMD / PERIMI16 ==================
// Генерируем набор диагоналей (off = j-i) вокруг band_centers и собираем всё в PerimD.
// Нулей нет: учитываем только те диагонали, у которых длина>0.
typedef struct {
    int off;
    int len;
    int pos; // позиция контакта на периметре
} DiagDesc;

static PerimD* perim_build_double(int n){
    // грубая верхняя оценка количества диагоналей
    int max_diags = band_count * (2*band_half_width+1);
    DiagDesc *diags = (DiagDesc*)malloc(sizeof(DiagDesc)*max_diags);
    int diag_cnt = 0;
    size_t total_len = 0;

    // собираем диагонали
    for(int bc=0; bc<band_count; ++bc){
        int c = band_centers[bc];
        for(int d=-band_half_width; d<=band_half_width; ++d){
            int off = c + d;
            if(off <= -n+1 || off >= n) continue; // слишком далеко, нет пересечения
            int abs_off = (off<0)?-off:off;
            int len = n - abs_off; // длина диагонали
            if(len<=0) continue;
            int pos = perimeter_contact_pos(n, off);
            diags[diag_cnt].off = off;
            diags[diag_cnt].len = len;
            diags[diag_cnt].pos = pos;
            total_len += (size_t)len;
            diag_cnt++;
        }
    }
    if(diag_cnt==0){
        fprintf(stderr,"No diagonals generated!\n");
        exit(1);
    }

    // сортировка диагоналей по pos (простая O(k^2), k << n)
    for(int i=0;i<diag_cnt;++i){
        for(int j=i+1;j<diag_cnt;++j){
            if(diags[j].pos < diags[i].pos){
                DiagDesc tmp=diags[i]; diags[i]=diags[j]; diags[j]=tmp;
            }
        }
    }

    // выделяем PerimD
    PerimD *P = (PerimD*)calloc(1,sizeof(PerimD));
    if(!P){ fprintf(stderr,"alloc PerimD failed\n"); exit(1); }
    P->n = n;
    P->strips = diag_cnt;
    P->meta = (int*)malloc(sizeof(int)*2*diag_cnt);
    P->data = (double*)malloc(sizeof(double)*total_len);
    if(!P->meta || !P->data){
        fprintf(stderr,"alloc PerimD fields failed\n"); exit(1);
    }

    // заполняем meta и data
    int prev_pos = 0;
    size_t w = 0;
    for(int s=0; s<diag_cnt; ++s){
        int off = diags[s].off;
        int len = diags[s].len;
        int pos = diags[s].pos;
        int delta = (s==0) ? pos : (pos - prev_pos);
        prev_pos = pos;

        P->meta[2*s+0] = delta;
        P->meta[2*s+1] = len;

        // для диагонали off: i идёт от i0 до i0+len-1
        int i0 = (off<0)? -off : 0;
        for(int t=0;t<len;++t){
            int i = i0 + t;
            int j = i + off;
            // генерируем значение A[i,j] в (-1,1)
            double r = ((double)rand()/(double)RAND_MAX)*2.0 - 1.0;
            // усилим главную диагональ
            if(j==i) r = 0.5*r + (r>=0 ? 0.5 : -0.5);
            if(r>=1.0) r=0.999999;
            if(r<=-1.0) r=-0.999999;
            P->data[w++] = r;
        }
    }

    free(diags);
    return P;
}

static PerimI16* perim_from_double(const PerimD *Pd, int S){
    PerimI16 *Pi = (PerimI16*)calloc(1,sizeof(PerimI16));
    if(!Pi){ fprintf(stderr,"alloc PerimI16 failed\n"); exit(1); }
    Pi->n = Pd->n;
    Pi->strips = Pd->strips;
    size_t meta_len = (size_t)Pd->strips*2;
    // считаем количество элементов data
    // (можно суммировать meta[2*s+1])
    size_t total_len=0;
    for(int s=0; s<Pd->strips; ++s) total_len += (size_t)Pd->meta[2*s+1];

    Pi->meta = (int*)malloc(sizeof(int)*meta_len);
    Pi->data = (int16_t*)malloc(sizeof(int16_t)*total_len);
    if(!Pi->meta || !Pi->data){
        fprintf(stderr,"alloc PerimI16 fields failed\n"); exit(1);
    }

    for(size_t i=0;i<meta_len;++i) Pi->meta[i]=Pd->meta[i];
    for(size_t i=0;i<total_len;++i){
        double scaled = Pd->data[i]*(double)S;
        Pi->data[i] = (int16_t)scaled; // предполагается, что не переполняется
    }
    return Pi;
}

// ================== ВАШИ ФУНКЦИИ: perim_matvec_* ==================
// Здесь — простые "референсы", которые можно заменить на свою реализацию.
// Формат входа: P->meta = [delta_up, width], P->data — значения диагоналей подряд.

void perim_matvec_double_ref(const PerimD *P, const double *x, double *y){
    int n = P->n;
    for(int i=0;i<n;++i) y[i]=0.0;

    int pos = 0;
    size_t rp = 0;
    for(int s=0; s<P->strips; ++s){
        int delta_up = P->meta[2*s+0];
        int width    = P->meta[2*s+1];
        pos += delta_up;

        // восстановим off из позиции на периметре
        int off = (pos < n) ? (pos - (n - 1)) : (pos - n);
        int i0  = (off<0)? -off : 0;

        for(int t=0; t<width; ++t){
            int i = i0 + t;
            int j = i + off;
            double a = P->data[rp++];
            y[i] += a * x[j];
        }
    }
}

void perim_matvec_i16_ref(const PerimI16 *P, const int16_t *xS, int32_t *yS2){
    int n = P->n;
    for(int i=0;i<n;++i) yS2[i]=0;

    int pos = 0;
    size_t rp = 0;
    for(int s=0; s<P->strips; ++s){
        int delta_up = P->meta[2*s+0];
        int width    = P->meta[2*s+1];
        pos += delta_up;

        int off = (pos < n) ? (pos - (n - 1)) : (pos - n);
        int i0  = (off<0)? -off : 0;

        for(int t=0; t<width; ++t){
            int i = i0 + t;
            int j = i + off;
            int32_t a = (int32_t)P->data[rp++];
            int32_t b = (int32_t)xS[j];
            yS2[i] += a * b;  // внутренний цикл считает "кусок" одного y[i]
        }
    }
}

// ================== MAIN ==================
int main(void){
    long long seed = (long long)time(NULL);
    printf("seed: %lld\n", seed);
    srand((unsigned)seed);

    // 1) Строим PerimD (double)
    clock_t t_buildD = clock();
    PerimD *Pd = perim_build_double(n);
    t_buildD = clock() - t_buildD;

    // 2) Строим PerimI16 (из PerimD)
    clock_t t_mat_quant = clock();
    PerimI16 *Pi = perim_from_double(Pd, S);
    t_mat_quant = clock() - t_mat_quant;

    // Векторы/буферы
    double  *x         = (double*) malloc(sizeof(double)*n);
    double  *y_d       = (double*) malloc(sizeof(double)*n);
    double  *y_q_perim = (double*) malloc(sizeof(double)*n);
    double  *tmp       = (double*) malloc(sizeof(double)*n);

    int16_t *xS        = (int16_t*)malloc(sizeof(int16_t)*n);
    int16_t *xS2       = (int16_t*)malloc(sizeof(int16_t)*n);
    int32_t *yS2       = (int32_t*)malloc(sizeof(int32_t)*n);

    if(!x || !y_d || !y_q_perim || !tmp || !xS || !xS2 || !yS2){
        fprintf(stderr,"alloc vectors failed\n");
        return 1;
    }

    // случайный x в (-1,1)
    for(int i=0;i<n;++i){
        double rx = ((double)rand()/(double)RAND_MAX)*2.0 - 1.0;
        if(rx>=1.0) rx=0.999999;
        if(rx<=-1.0) rx=-0.999999;
        x[i]=rx;
    }

    // ----- baseline: perim_matvec_double_ref -----
    clock_t t_double = 0;
    for(int e=0; e<epochs; ++e){
        for(int i=0;i<n;++i) tmp[i]=x[i];
        clock_t t0=clock();
        for(int it=0; it<iterationsCount; ++it){
            perim_matvec_double_ref(Pd, tmp, y_d);
            for(int i=0;i<n;++i) tmp[i]=y_d[i];
        }
        t_double += clock()-t0;
    }

    // ----- optimized: perim_matvec_i16_ref -----
    clock_t t_conv     = 0;
    clock_t t_matmul_i = 0;

    for(int e=0; e<epochs; ++e){
        clock_t t0=clock();
        quantize_x_i16(x, n, xS, S);
        t_conv += clock()-t0;

        for(int it=0; it<iterationsCount; ++it){
            t0=clock();
            perim_matvec_i16_ref(Pi, xS, yS2);
            t_matmul_i += clock()-t0;

            t0=clock();
            dequantize_y_update_xS(yS2, y_q_perim, xS2, n, S);
            for(int i=0;i<n;++i) xS[i]=xS2[i];
            t_conv += clock()-t0;
        }
    }

    // ----- ошибка на одной итерации -----
    perim_matvec_double_ref(Pd, x, y_d);
    quantize_x_i16(x, n, xS, S);
    perim_matvec_i16_ref(Pi, xS, yS2);
    dequantize_y_update_xS(yS2, y_q_perim, xS2, n, S);

    double err_l1 = rel_l1(y_d, y_q_perim, n);
    double err_l2 = rel_l2(y_d, y_q_perim, n);

    // ----- вывод -----
    double tD   = (double)t_double   / CLOCKS_PER_SEC / epochs;
    double tMM  = (double)t_matmul_i / CLOCKS_PER_SEC / epochs;
    double tCV  = (double)t_conv     / CLOCKS_PER_SEC / epochs;
    double tBD  = (double)t_buildD   / CLOCKS_PER_SEC;
    double tAQ  = (double)t_mat_quant/ CLOCKS_PER_SEC;

    printf("n = %d, strips = %d\n", n, Pd->strips);
    printf("S = %d\n", S);

    printf("RelError (perim-int16) L1: %.6f%%\n", err_l1);
    printf("RelError (perim-int16) L2: %.6f%%\n", err_l2);

    printf("Time (perim-double baseline):          %.6f s/epoch\n", tD);
    printf("Time (perim-int16 matmul only):        %.6f s/epoch\n", tMM);
    printf("Time (perim-int16 quant/dequant epoch):%.6f s/epoch\n", tCV);
    printf("Time (perim-double build once):        %.6f s\n",      tBD);
    printf("Time (perim-int16  build once):        %.6f s\n",      tAQ);

    if(tMM>0.0) printf("Speedup (double / perim-int16):       %.2fx\n", tD/tMM);

    // суммы результатов
    double sum_y_d=0.0, sum_y_q=0.0;
    for(int i=0;i<n;++i){
        sum_y_d += y_d[i];
        sum_y_q += y_q_perim[i];
    }
    printf("Sum(y_perim double):       %.10e\n", sum_y_d);
    printf("Sum(y_perim int16 dequant):%.10e\n", sum_y_q);

    // очистка
    free(x); free(y_d); free(y_q_perim); free(tmp);
    free(xS); free(xS2); free(yS2);
    perimd_free(Pd);
    perimi16_free(Pi);
    return 0;
}
