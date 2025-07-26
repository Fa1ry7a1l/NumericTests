#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define max(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a > _b ? _a : _b; })

#define min(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a < _b ? _a : _b; })

#include <time.h>

int M = 16;

double *a;
unsigned short *a_int;
double *b;
unsigned short *b_int;
unsigned int *c_int;

double *res;
double *res2;


long long time_spent_usual = 0;
long long time_spent_optimized = 0;


// void generateData() {
//     for (int i = 0; i < n1; i++) {
//         a[i] = rand() / (double) ((long) RAND_MAX + 1);
//         a_int[i] = a[i] * M;
//     }
//     for (int i = 0; i < n2; i++) {
//         b[i] = rand() / (double) ((long) RAND_MAX + 1);
//         b_int[i] = b[i] * M;
//     }
//     for (int i = 0; i < n1 + n2; i++) {
//         c_int[i] = 0;
//         res2[i] = 0;
//         res[i] = 0;
//     }
//
//
//     time_spent_optimized = 0;
//     time_spent_usual = 0;
// }
//
// void print_errors(void) {
//     double delta1 = 0, delta2 = 0, norm1 = 0, norm2 = 0, relativeError1 = 0, relativeError2 = 0;
//
//     for (int i = 0; i < n1 + n2; i++) {
//         delta1 += fabs(res[i] - res2[i]);
//         norm1 += fabs(res[i]);
//
//         delta2 += pow(res[i] - res2[i], 2);
//         norm2 += pow(res[i], 2);
//     }
//
//     relativeError1 = delta1 * (double) 100 / norm1;
//     relativeError2 = sqrt(delta2) * (double) 100 / sqrt(norm2);
//
//     printf("Relative Error1: %f norm: %f\n", relativeError1, norm1);
//     printf("Relative Error2: %f norm: %f\n", relativeError2, norm2);
// }

void banded_matvec_mult_int(long *Ab, long *x, long *y, int n, int kl, int ku, long *tmp, int n0) {
    for (int i = 0; i < n; ++i) {
        y[i] = 0;

    //     // j: столбцы от i - kl до i + ku
    int j_start = (i - kl) > 0 ? (i - kl) : 0;
    int j_end = (i + ku) < (n - 1) ? (i + ku) : (n - 1);
    int m = (j_end - j_start + 1) / n0;
    int m_up = (j_end - j_start + 1 + n0 - 1) / n0;
    for (int j = 0; j < m_up; j++) {
        tmp[j] = 0;
    }
    //основная часть
    for (int k1 = 0; k1 < m; k1++) {
        int t = k1 * n0;
        for (int k2 = 0; k2 < n0; k2++) {
            int j = t + k2 + j_start;
            int band_row = ku + i - j;
            if (j >= n || j < 0) {
                printf("error in main i=%d m=%d, k1=%d, k2=%d\n, t=%d, j_start=%d", i, m, k1, k2, t, j_start);
            }
            tmp[k1] += Ab[band_row * n + j] * x[j];
        }
    }

    int k1_2 = m;
    int t = m * n0;
    for (int k2 = 0; k2 < (j_end - j_start + 1) % n0; k2++) {
        int j = t + k2 + j_start;
        int band_row = ku + i - j;
        if (j >= n || j < 0) {
            printf("error in mod i=%d m=%d, k2=%d\n, t=%d, j_start=%d", i, m, k2, t, j_start);
        }
        tmp[k1_2] += Ab[band_row * n + j] * x[j];
    }

    for (int j = 0; j < m_up; j++) {
        y[i] += tmp[j];
    }
    }
}

// Умножение ленточной матрицы A на вектор x: y = A * x
// A хранится как одномерный массив в ленточном формате
// n — размер матрицы, k1 — нижняя полуширина, k2 — верхняя полуширина
void banded_matvec_mult(double *Ab, double *x, double *y, int n, int kl, int ku) {
    for (int i = 0; i < n; ++i) {
        y[i] = 0.0;

        // j: столбцы от i - kl до i + ku
        int j_start = (i - kl) > 0 ? (i - kl) : 0;
        int j_end = (i + ku) < (n - 1) ? (i + ku) : (n - 1);

        for (int j = j_start; j <= j_end; ++j) {
            int band_row = ku + i - j;
            y[i] += Ab[band_row * n + j] * x[j];
        }
    }
}


double randd() {
    return rand() / (double) ((long) RAND_MAX + 1);
}

int main() {
    int n = 200; // размерность
    int k1 = 100; // нижняя полуширина
    int k2 = 100; // верхняя полуширина
    int n0 = 20;
    int m = k1 - k2 + 1;

    // Выделим память
    double *A = (double *) malloc(n * (k1 + k2 + 1) * sizeof(double));
    long *Ashort = (long *) malloc(n * (k1 + k2 + 1) * sizeof(long));
    double *x = (double *) malloc(n * sizeof(double));
    long *xshort = (long *) malloc(n * sizeof(long));
    double *y = (double *) malloc(n * sizeof(double));
    long *y2 = (long *) malloc(n * sizeof(long));
    long *tmp = (long *) malloc(n * sizeof(long));

    long long seed = time(NULL);
    printf("seed: %lld\n", seed);

    srand(seed);

    for (int i = 0; i < m + 1; i++) {
        tmp[i] = 0;
    }

    // Заполним матрицу и вектор
    for (int i = 0; i < n; i++) {
        x[i] = randd();
        xshort[i] = x[i] * M;
        y[i] = 0;
        y2[i] = 0;
        for (int j = -k1; j <= k2; j++) {
            int col = i + j;
            if (col >= 0 && col < n) {
                A[i * (k1 + k2 + 1) + (j + k1)] =  randd();// Пример: все элементы в ленте — 1.0
                Ashort[i * (k1 + k2 + 1) + (j + k1)] = A[i * (k1 + k2 + 1) + (j + k1)]*M;
            } else {
                A[i * (k1 + k2 + 1) + (j + k1)] = 0.0;
                Ashort[i * (k1 + k2 + 1) + (j + k1)] = 0;
            }
        }
    }



    // Умножение
    banded_matvec_mult(A, x, y, n, k1, k2);
    banded_matvec_mult_int(Ashort, xshort, y2, n, k1, k2,tmp,n0);

    // Вывод результата
    printf("y = [");
    for (int i = 0; i < n; i++) {
        printf(" %.2f", y[i]);
    }
    printf(" ]\n");

    printf("y2 = [");
    for (int i = 0; i < n; i++) {
        printf(" %.2f", y2[i] / ((double) M * M));
    }
    printf(" ]\n");


    // Освобождение памяти
    free(y);
    free(x);
    free(A);
    free(y2);
    free(tmp);

    free(Ashort);
    free(xshort);


    return 0;
}
