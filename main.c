#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// ===== ПАРАМЕТРЫ ЗАДАЁМ РУКАМИ ЗДЕСЬ =====
const int n  = 127000;   // длина векторов
const int n0 = 128;    // размер блока
const int M  = 2048;  // масштаб для квантования в short
// =========================================

const int m = (n+n0-1)/n0; // число блоков = n / n0 (может быть хвост)

double *a, *b;
float *af, *bf;
unsigned short *a1, *b1;
unsigned long long *temp;

double res_double, res_short_block, res_short_long;
float  res_float;

float epoh = 20.0f;

// ======== DOUBLE: скалярное произведение ========
double calculate_double_time() {
    clock_t begin = clock();
    double s = 0.0;
    for (int i = 0; i < n; i++)
        s += a[i] * b[i];
    res_double = s;
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

// ======== FLOAT: скалярное произведение ========
double calculate_float_time() {
    clock_t begin = clock();
    float s = 0.0f;
    for (int i = 0; i < n; i++)
        s += af[i] * bf[i];
    res_float = s;
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

// ======== SHORT (блочная логика) ========
double calculate_short_block_time() {
    clock_t begin = clock();
    unsigned long long S = 0;

    // блоки длиной n0
    for (int i = 0; i < m; i++) {
        int t = i * n0;
        unsigned long long local = 0;
        for (int j = 0; j < n0 && t < n; j++) {
            local += (unsigned long long)a1[t] * (unsigned long long)b1[t];
            t++;
        }
        temp[i] = local;
    }

    for (int i = 0; i < m; i++)
        S += temp[i];

    // хвост, если n не кратно n0
    int start_tail = m * n0;
    for (int i = start_tail; i < n; i++)
        S += (unsigned long long)a1[i] * (unsigned long long)b1[i];

    res_short_block = (double)S / ((double)M * (double)M);
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

// ======== SHORT (одним большим проходом) ========
double calculate_short_long_time() {
    clock_t begin = clock();
    unsigned long long S = 0;

    for (int i = 0; i < n; i++)
        S += (unsigned long long)a1[i] * (unsigned long long)b1[i];

    res_short_long = (double)S / ((double)M * (double)M);
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

// ======== генерация данных ========
void generateData() {
    a  = malloc(sizeof(double) * n);
    b  = malloc(sizeof(double) * n);
    af = malloc(sizeof(float)  * n);
    bf = malloc(sizeof(float)  * n);
    a1 = malloc(sizeof(unsigned short) * n);
    b1 = malloc(sizeof(unsigned short) * n);
    temp = malloc(sizeof(unsigned long long) * m);

    if (!a || !b || !af || !bf || !a1 || !b1 || !temp) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        double x = rand() / (double)(RAND_MAX + 1);
        double y = rand() / (double)(RAND_MAX + 1);

        a[i] = x;
        b[i] = y;

        af[i] = (float)x;
        bf[i] = (float)y;

        unsigned int qx = (unsigned int)(x * M);
        unsigned int qy = (unsigned int)(y * M);
        if (qx >= (unsigned int)M) qx = M - 1;
        if (qy >= (unsigned int)M) qy = M - 1;

        a1[i] = (unsigned short)qx;
        b1[i] = (unsigned short)qy;
    }
}

int main(void) {
    printf("%d\n",sizeof(float));

    printf("%d\n",sizeof(double));
    srand((unsigned int)time(NULL));

    generateData();

    double Tdouble = 0.0, Tfloat = 0.0, TshortB = 0.0, TshortL = 0.0;
    double errF_abs = 0.0, errF_rel = 0.0;
    double errB_abs = 0.0, errB_rel = 0.0;
    double errL_abs = 0.0, errL_rel = 0.0;

    // ======= ЧАСТЬ 1 — ЭПОХИ ТОЛЬКО ДЛЯ double =======
    for (int z = 0; z < (int)epoh; z++) {
        Tdouble += calculate_double_time();
    }

    // ======= ЧАСТЬ 2 — ЭПОХИ ТОЛЬКО ДЛЯ float =======
    for (int z = 0; z < (int)epoh; z++) {

        Tfloat += calculate_float_time();

        // Ошибка float только относительно double
        double denom = fabs(res_double);
        double df = (double)res_float - res_double;

        errF_abs += fabs(df);
        errF_rel += fabs(df) / denom;
    }

    // ======= ЧАСТЬ 3 — ЭПОХИ ТОЛЬКО ДЛЯ short =======
    for (int z = 0; z < (int)epoh; z++) {

        TshortB += calculate_short_block_time();
        TshortL += calculate_short_long_time();

        double denom = fabs(res_double) + 1e-30;

        double db = res_short_block - res_double;
        double dl = res_short_long  - res_double;

        errB_abs += fabs(db);
        errB_rel += fabs(db) / denom;

        errL_abs += fabs(dl);
        errL_rel += fabs(dl) / denom;
    }

    // усреднение по эпохам
    Tdouble /= epoh;
    Tfloat  /= epoh;
    TshortB /= epoh;
    TshortL /= epoh;

    errF_abs /= epoh; errF_rel /= epoh;
    errB_abs /= epoh; errB_rel /= epoh;
    errL_abs /= epoh; errL_rel /= epoh;

    printf("n=%d; n0=%d; M=%d\n", n, n0, M);
    printf("timeDouble=%f s\n",      Tdouble);
    printf("timeFloat=%f s\n",       Tfloat);
    printf("timeShortBlock=%f s\n",  TshortB);
    printf("timeShortLong=%f s\n\n", TshortL);

    printf("errFloatAbs=%e;  errFloatRel=%e\n",       errF_abs, errF_rel);
    printf("errShortBlockAbs=%e; errShortBlockRel=%e\n", errB_abs, errB_rel);
    printf("errShortLongAbs=%e;  errShortLongRel=%e\n",  errL_abs, errL_rel);

    free(a); free(b);
    free(af); free(bf);
    free(a1); free(b1);
    free(temp);

    return 0;
}
