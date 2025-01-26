#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

int n = 1000;
int n0 = 255;
int M = 4096;
int M1 = 4096;
int M2 = 4096;
double **a;
double *b;
double *x;
double *res;
double *res2;
int iterationsCount = 100;
int epochs = 10;


long long time_spent_usual = 0;
long long time_spent_optimized = 0;

unsigned int *temp;
unsigned short *x_int;
unsigned short **a_int;
unsigned int *b_int;
unsigned long long *res_int;

double *x_float;


void generateMatrix(double **a) {
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) {
            a[i][j] = rand() / (double) ((long) RAND_MAX + 1);
            sum += a[i][j];
        }
        sum *= 1.36;

        for (int j = 0; j < n; j++) //todo не должно быть нормальзации
        {
            a[i][j] /= sum;
        }
    }
}

void free_ram(void) {
    for (int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(b);
    free(res);
    free(x);
}

void calculate_floats(void) {
    for (int i = 0; i < n; i++) {
        x_float[i] = x[i];
    }

    time_t begin = clock();

    for (int it = 0; it < iterationsCount; it++) {
        for (int i = 0; i < n; i++) {
            res[i] = 0;
            for (int j = 0; j < n; j++) {
                res[i] += a[i][j] * x_float[j];
            }
            res[i] += b[i];
        }

        for (int i = 0; i < n; i++) {
            x_float[i] = res[i];
        }
    }
    time_spent_usual = clock() - begin;

    //free(x_float);
}

void calculate_ints(void) {
    for (int i = 0; i < n; i++) //домножаем вектор x на M
    {
        x_int[i] = x[i] * M;
    }
    for (int i = 0; i < n; i++) //домножаем матрицу а на M
    {
        a_int[i] = malloc(sizeof(unsigned short) * n);
        for (int j = 0; j < n; j++) {
            a_int[i][j] = a[i][j] * M;
        }
    }
    for (int i = 0; i < n; i++) //домножаем вектор b на M*М
    {
        b_int[i] = b[i] * M * M;
    }
    for (int i = 0; i < n; i++) {
        res_int[i] = 0;
    }


    time_t begin = clock();

    for (int it = 0; it < iterationsCount; it++) {
        for (int i = 0; i < n; i++) {
            unsigned long long sum = 0;

            const int m = n / n0;
            for (int j = 0; j < m; j++) {
                temp[j] = 0;
                int t = j * n0;
                for (int k = 0; k < n0; k++, t++) {
                    temp[j] += a_int[i][t] * x_int[t];
                }
            }

            for (int j = 0; j < m; j++) {
                sum += temp[j];
            }
            temp[0] = 0;
            for (int j = n0 * m; j < n; j++) {
                temp[0] += a_int[i][j] * x_int[j];
            }
            sum += temp[0];


            sum += b_int[i];
            //посчитали в целых числах перемножение строки матрицы на вектор x + соответствующий элемент b
            res2[i] = sum / (double) M; //делим на M
        }

        for (int i = 0; i < n; i++) {
            x_int[i] = res2[i];
            //переприсываиваем результаты вычисления во входной вектор, чтобы пойти на новую итерацию.
        }
    }

    time_spent_optimized = clock() - begin;

    for (int i = 0; i < n; i++) {
        res2[i] = res2[i] / (double) M; //делим результирующий вектор еще раз на M
    }


    // free(x_int);
    // // for (int i = 0; i < n; i++) {
    // //     free(a_int[i]);
    // // }
    // free(a_int);
    // free(b_int);
    // free(res_int);
    // free(temp);
}

void print_errors(void) {
    double delta1 = 0, delta2 = 0, norm1 = 0, norm2 = 0, relativeError1 = 0, relativeError2 = 0;

    for (int i = 0; i < n; i++) {
        delta1 += fabs(res[i] - res2[i]);
        norm1 += fabs(res[i]);

        delta2 += pow(res[i] - res2[i], 2);
        norm2 += pow(res[i], 2);
    }

    relativeError1 = delta1 * (double) 100 / norm1;
    relativeError2 = sqrt(delta2) * (double) 100 / sqrt(norm2);

    printf("Relative Error1: %f norm: %f\n", relativeError1, norm1);
    printf("Relative Error2: %f norm: %f\n", relativeError2, norm2);
}

void check_difference() {
    for (int i = 0; i < n; i++) {
        if (res[i] != res2[i]) {
            printf("%f %f %f\n", res[i], res2[i], res[i] - res2[i]);
        }
    }
}

int main(void) {
    temp = (unsigned int *) malloc(sizeof(unsigned int) * n0);
    x_int = (unsigned short *) malloc(sizeof(unsigned short) * n);
    a_int = (unsigned short **) malloc(sizeof(unsigned short *) * n);
    b_int = (unsigned int *) calloc(n, sizeof(unsigned int));
    res_int = (unsigned long long *) malloc(sizeof(unsigned long long) * n);
    x_float = (double *) malloc(sizeof(double) * n);
    long long seed = time(NULL);
    printf("seed: %lld\n", seed);
    srand(seed);

    //for (int z = 128;z < 5000;z*=2) {
    //    n=z;


    a = (double **) malloc(sizeof(double *) * n);
    for (int i = 0; i < n; i++) {
        a[i] = (double *) malloc(sizeof(double) * n);
        for (int j = 0; j < n; j++) {
            a[i][j] = 0;
        }
    }
    generateMatrix(a);
    b = (double *) malloc(sizeof(double) * n);
    x = (double *) malloc(sizeof(double) * n);
    res = (double *) malloc(sizeof(double) * n);
    res2 = (double *) malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        b[i] = rand() / (double) ((long) RAND_MAX + 1);
        x[i] = rand() / (double) ((long) RAND_MAX + 1);
        res[i] = 0;
        res2[i] = 0;
    }

    for (int i = 0; i < epochs; i++) {
        calculate_floats();
        calculate_ints();
        if (i % 10 == 0)
            printf("IterationsCount: %d\n", i);
        for (int j = 0; j < n; j++) {
            res[i] = res2[i] = 0;
        }
    }


    print_errors();
    //check_difference();

    printf("n=%d\n", n);
    printf("Time spent optimized: %f\n", (double) time_spent_optimized / CLOCKS_PER_SEC / epochs);
    printf("Time spent usual: %f\n", (double) time_spent_usual / CLOCKS_PER_SEC / epochs);

    //free_ram();
    //}
    return 0;
}
