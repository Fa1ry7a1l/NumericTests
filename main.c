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

int n1 = 32768;
int n2 = 32768;
int n0 = 40;
int M = 1000;

double *a;
unsigned short *a_int;
double *b;
unsigned short *b_int;
unsigned int *c_int;

double *res;
double *res2;


long long time_spent_usual = 0;
long long time_spent_optimized = 0;


void calculate_floats(void) {
    time_t begin = clock();
    for (int i = 0; i < n1 + n2-1; i++) {
        int L = max(0, i-n2 + 1);
        int R = min(i, n1 - 1);
        for (int j = L; j <= R; j++) {
            res[i] += a[j] * b[i - j];
        }
    }
    time_spent_usual = (clock() - begin);
}

void sum(int n, int res, int *a) {
    for (int i = 0; i < n; i++) {
        res += a[i];
    }
}

//full sum optimization
void sum_opt(int n, int n0, int res, int *temp_res, int *a) {
    int m = n / n0;
    for (int i = 0; i < m; i++) {
        temp_res[i] = 0;
        int t = i * n0;
        for (int j = 0; j < n0; j++) {
            temp_res[i] += a[j + t];
        }
    }
    for (int i = 0; i < m; i++) {
        res += temp_res[i];
    }

    temp_res[0] = 0;
    int round_n = m * n0;
    for (int i = round_n; i < n; i++) {
        temp_res[0] += a[i];
    }
    res += temp_res[0];
}

void calculate_ints(void) {
    time_t begin = clock();
    for (int i = 0; i < n1 + n2-1; i++) {
        int L = max(0, i-n2 + 1);
        int R = min(i, n1 - 1);

        //part for sum
        int n = R - L + 1;
        int m = n / n0;
        for (int j = 0; j < m; j++) {
            c_int[j] = 0;
            int t = j * n0 + L;
            for (int k = 0; k < n0; k++) {
                c_int[j] += a_int[k + t] * b_int[i - (k + t)];
            }
        }
        for (int j = 0; j < m; j++) {
            res2[i] += c_int[j];
        }
        c_int[0] = 0;
        for (int j = m * n0; j < n; j++) {
            c_int[0] += a_int[j] * b_int[i - j];
        }
        res2[i] += c_int[0];
        res2[i] /= M * M;
    }


    time_spent_optimized = (clock() - begin);
}

void generateData() {
    for (int i = 0; i < n1; i++) {
        a[i] = rand() / (double) ((long) RAND_MAX + 1);
        a_int[i] = a[i] * M;
    }
    for (int i = 0; i < n2; i++) {
        b[i] = rand() / (double) ((long) RAND_MAX + 1);
        b_int[i] = b[i] * M;
    }
    for (int i = 0; i < n1 + n2; i++) {
        c_int[i] = 0;
        res2[i] = 0;
        res[i] = 0;
    }


    time_spent_optimized = 0;
    time_spent_usual = 0;
}

void print_errors(void) {
    double delta1 = 0, delta2 = 0, norm1 = 0, norm2 = 0, relativeError1 = 0, relativeError2 = 0;

    for (int i = 0; i < n1 + n2; i++) {
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

int main(void) {
    a = (double *) malloc(sizeof(double) * n1);
    a_int = (unsigned short *) malloc(sizeof(unsigned short) * n1);
    b = (double *) malloc(sizeof(double) * n2);
    b_int = (unsigned short *) malloc(sizeof(unsigned short) * n2);
    c_int = (unsigned int *) malloc(sizeof(unsigned int) * (n1 + n2));
    res = (double *) malloc(sizeof(double) * (n1 + n2));
    res2 = (double *) malloc(sizeof(double) * (n1 + n2));

    long long seed = time(NULL);
    printf("seed: %lld\n", seed);

    srand(seed);

    printf("n1=%d, n0=%d, M=%d\n", n1, n0, M);


    generateData();
    calculate_floats();
    calculate_ints();


    print_errors();

    printf("n1=%d n2=%d  m=%d\n", n1, n2, (n1 + n2) / n0);
    printf("for usual time = %lld\n", time_spent_usual);
    printf("for optimized time = %lld\n", time_spent_optimized);
    printf("\n");


    free(a);
    free(a_int);
    free(b);
    free(b_int);
    free(c_int);
    free(res);
    free(res2);


    return 0;
}
