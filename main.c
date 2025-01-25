#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

int n = 10000;
int *a;
int *b;
int *res;
int *res2;
int epochs = 10;


long long time_spent_usual = 0;
long long time_spent_optimized = 0;


void free_ram(void) {
    free(a);
    free(b);
    free(res);
    free(res2);
}

void calculate_floats(void) {
    time_t begin = clock();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i] += a[i] * b[(j - i + n) % n];
        }
    }
    time_spent_usual = clock() - begin;
}

void calculate_ints(void) {
    time_t begin = clock();

    for (int i = 0; i < n; i++) {
        int j;
        for (j = 0; j < i; j++) {
            res2[i] += a[i] * b[j - i + n];
        }
        for (j = i; j < n; j++) {
            res2[i] += a[i] * b[j - i];
        }
    }
    time_spent_optimized = clock() - begin;
}


int main(void) {
    long long seed = time(NULL);
    printf("seed: %lld\n", seed);
    srand(seed);


    a = (int *) malloc(sizeof(int *) * n);
    b = (int *) malloc(sizeof(int) * n);
    res = (int *) malloc(sizeof(int) * n);
    res2 = (int *) malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) {
        b[i] = rand() /*/ (double) ((long) RAND_MAX + 1)*/;
        a[i] = rand() /*/ (double) ((long) RAND_MAX + 1)*/;
        res[i] = 0;
        res2[i] = 0;
    }

    for (int i = 0; i < epochs; i++) {
        calculate_floats();
        calculate_ints();
        if (i % 10 == 0)
            printf("IterationsCount: %d\n", i);
    }


    printf("n=%d\n", n);
    printf("Time spent optimized: %f\n", (double) time_spent_optimized / CLOCKS_PER_SEC / epochs);
    printf("Time spent usual: %f\n", (double) time_spent_usual / CLOCKS_PER_SEC / epochs);

    free_ram();

    return 0;
}
