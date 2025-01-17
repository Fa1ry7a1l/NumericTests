#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

int n = 128;
int n0 = 128;
int m = 100;
int M = 2048;
double *a;
int *a1;
int *temp;

double res;
double res2;


long long time_spent_usual = 0;
long long time_spent_optimized = 0;


void calculate_floats(void) {
    time_t begin = clock();
    for (int i = 0; i < n; i++) {
        res += a[i];
    }
    time_spent_usual = (clock() - begin) ;
}

void calculate_ints(void) {
    time_t begin = clock();
    long unsigned long S = 0;
    for (int i = 0; i < m; i++) {
        int t = i * n0;
        for (int j = 0; j < n0; j++) {
            temp[i] += a1[t++];
        }
    }
    for (int i = 0; i < m; i++) {
        S += temp[i];
    }
    res2 = (double) S / (double) M;
    time_spent_optimized = (clock() - begin) ;
}

void generateData() {
    a = (double *) malloc(sizeof(double) * n);
    a1 = (int *) malloc(sizeof(int) * n);
    temp = (int *) malloc(sizeof(int) * m);

    for (int i = 0; i < m; i++) {
        temp[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        a[i] = rand() / (double) ((long) RAND_MAX + 1);
        a1[i] = a[i] * M;
    }

    res = 0;
    res2 = 0;
    time_spent_optimized = 0;
    time_spent_usual = 0;
}

int main(void) {
    long long seed = time(NULL);
    printf("seed: %lld\n", seed);

    srand(seed);

    for (int i = 100; i < 1000000000; i *= 10) {
        n = i;
        for (int j = 64; j < n; j *= 2) {
            M = j;
            for (int k = 1; k <= n / M; k *= 2) {
                n0 = k;
                m = n / n0;
                generateData(n);
                calculate_floats();
                calculate_ints();

                free(a);
                free(a1);
                free(temp);
                printf("n=%d n0=%d  M=%d\n", n, n0, M);
                printf("for usual time = %lld\n",time_spent_usual);
                printf("for optimized time = %lld\n",time_spent_optimized);
                printf("\n");
            }



        }
    }


    return 0;
}
