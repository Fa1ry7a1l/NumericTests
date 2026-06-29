#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static size_t n = 127000;
static size_t n0 = 128;
static const int M = 2048;
static size_t m = 0;
static int epochs = 20;

static double *a, *b;
static float *af, *bf;
static unsigned short *a1, *b1;
static unsigned long long *temp;

static double res_double, res_short_block, res_short_long;
static float res_float;

static double bytes_to_gib(size_t bytes) {
    return (double)bytes / (1024.0 * 1024.0 * 1024.0);
}

static size_t checked_add(size_t left, size_t right, const char *label) {
    if (SIZE_MAX - left < right) {
        fprintf(stderr, "size overflow while estimating %s\n", label);
        exit(2);
    }
    return left + right;
}

static size_t checked_mul(size_t count, size_t item_size, const char *label) {
    if (item_size != 0 && count > SIZE_MAX / item_size) {
        fprintf(stderr, "size overflow while estimating %s\n", label);
        exit(2);
    }
    return count * item_size;
}

static void *malloc_array(size_t count, size_t item_size, const char *label) {
    size_t bytes = checked_mul(count, item_size, label);
    void *ptr = malloc(bytes);
    if (!ptr) {
        fprintf(stderr, "malloc error for %s: %zu bytes (%.3f GiB)\n",
                label, bytes, bytes_to_gib(bytes));
        exit(1);
    }
    return ptr;
}

static size_t parse_size_arg(const char *text, const char *label) {
    char *end = NULL;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value == 0 ||
        value > (unsigned long long)SIZE_MAX) {
        fprintf(stderr, "bad %s value: %s\n", label, text);
        exit(2);
    }
    return (size_t)value;
}

static int parse_int_arg(const char *text, const char *label) {
    char *end = NULL;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value <= 0 ||
        value > INT_MAX) {
        fprintf(stderr, "bad %s value: %s\n", label, text);
        exit(2);
    }
    return (int)value;
}

static void configure(int argc, char **argv) {
    if (argc > 4) {
        fprintf(stderr, "usage: %s [n] [epochs] [n0]\n", argv[0]);
        exit(2);
    }

    if (argc > 1) {
        n = parse_size_arg(argv[1], "n");
    }
    if (argc > 2) {
        epochs = parse_int_arg(argv[2], "epochs");
    }
    if (argc > 3) {
        n0 = parse_size_arg(argv[3], "n0");
    }

    m = n / n0 + (n % n0 != 0);
}

static size_t estimated_memory_bytes(void) {
    size_t total = 0;
    total = checked_add(total, checked_mul(n, sizeof(*a), "a"), "total");
    total = checked_add(total, checked_mul(n, sizeof(*b), "b"), "total");
    total = checked_add(total, checked_mul(n, sizeof(*af), "af"), "total");
    total = checked_add(total, checked_mul(n, sizeof(*bf), "bf"), "total");
    total = checked_add(total, checked_mul(n, sizeof(*a1), "a1"), "total");
    total = checked_add(total, checked_mul(n, sizeof(*b1), "b1"), "total");
    total = checked_add(total, checked_mul(m, sizeof(*temp), "temp"), "total");
    return total;
}

static double calculate_double_time(void) {
    clock_t begin = clock();
    double s = 0.0;
    for (size_t i = 0; i < n; i++) {
        s += a[i] * b[i];
    }
    res_double = s;
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

static double calculate_float_time(void) {
    clock_t begin = clock();
    float s = 0.0f;
    for (size_t i = 0; i < n; i++) {
        s += af[i] * bf[i];
    }
    res_float = s;
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

static double calculate_short_block_time(void) {
    clock_t begin = clock();
    unsigned long long S = 0;

    for (size_t i = 0; i < m; i++) {
        size_t t = i * n0;
        unsigned long long local = 0;
        for (size_t j = 0; j < n0 && t < n; j++, t++) {
            local += (unsigned long long)a1[t] * (unsigned long long)b1[t];
        }
        temp[i] = local;
    }

    for (size_t i = 0; i < m; i++) {
        S += temp[i];
    }

    res_short_block = (double)S / ((double)M * (double)M);
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

static double calculate_short_long_time(void) {
    clock_t begin = clock();
    unsigned long long S = 0;

    for (size_t i = 0; i < n; i++) {
        S += (unsigned long long)a1[i] * (unsigned long long)b1[i];
    }

    res_short_long = (double)S / ((double)M * (double)M);
    return (double)(clock() - begin) / CLOCKS_PER_SEC;
}

static void generateData(void) {
    a = malloc_array(n, sizeof(*a), "a");
    b = malloc_array(n, sizeof(*b), "b");
    af = malloc_array(n, sizeof(*af), "af");
    bf = malloc_array(n, sizeof(*bf), "bf");
    a1 = malloc_array(n, sizeof(*a1), "a1");
    b1 = malloc_array(n, sizeof(*b1), "b1");
    temp = malloc_array(m, sizeof(*temp), "temp");

    for (size_t i = 0; i < n; i++) {
        double x = rand() / ((double)RAND_MAX + 1.0);
        double y = rand() / ((double)RAND_MAX + 1.0);

        a[i] = x;
        b[i] = y;

        af[i] = (float)x;
        bf[i] = (float)y;

        unsigned int qx = (unsigned int)(x * M);
        unsigned int qy = (unsigned int)(y * M);
        if (qx >= (unsigned int)M) {
            qx = M - 1;
        }
        if (qy >= (unsigned int)M) {
            qy = M - 1;
        }

        a1[i] = (unsigned short)qx;
        b1[i] = (unsigned short)qy;
    }
}

int main(int argc, char **argv) {
    configure(argc, argv);

    size_t estimated = estimated_memory_bytes();
    printf("sizeof(float)=%zu\n", sizeof(float));
    printf("sizeof(double)=%zu\n", sizeof(double));
    printf("n=%zu; n0=%zu; M=%d; epochs=%d\n", n, n0, M, epochs);
    printf("estimatedMemory=%zu bytes (%.3f GiB)\n", estimated,
           bytes_to_gib(estimated));

    srand((unsigned int)time(NULL));
    generateData();

    double Tdouble = 0.0, Tfloat = 0.0, TshortB = 0.0, TshortL = 0.0;
    double errF_abs = 0.0, errF_rel = 0.0;
    double errB_abs = 0.0, errB_rel = 0.0;
    double errL_abs = 0.0, errL_rel = 0.0;

    for (int z = 0; z < epochs; z++) {
        Tdouble += calculate_double_time();
    }

    for (int z = 0; z < epochs; z++) {
        Tfloat += calculate_float_time();

        double denom = fabs(res_double) + 1e-30;
        double df = (double)res_float - res_double;

        errF_abs += fabs(df);
        errF_rel += fabs(df) / denom;
    }

    for (int z = 0; z < epochs; z++) {
        TshortB += calculate_short_block_time();
        TshortL += calculate_short_long_time();

        double denom = fabs(res_double) + 1e-30;

        double db = res_short_block - res_double;
        double dl = res_short_long - res_double;

        errB_abs += fabs(db);
        errB_rel += fabs(db) / denom;

        errL_abs += fabs(dl);
        errL_rel += fabs(dl) / denom;
    }

    Tdouble /= epochs;
    Tfloat /= epochs;
    TshortB /= epochs;
    TshortL /= epochs;

    errF_abs /= epochs;
    errF_rel /= epochs;
    errB_abs /= epochs;
    errB_rel /= epochs;
    errL_abs /= epochs;
    errL_rel /= epochs;

    printf("timeDouble=%f s\n", Tdouble);
    printf("timeFloat=%f s\n", Tfloat);
    printf("timeShortBlock=%f s\n", TshortB);
    printf("timeShortLong=%f s\n\n", TshortL);

    printf("errFloatAbs=%e; errFloatRel=%e\n", errF_abs, errF_rel);
    printf("errShortBlockAbs=%e; errShortBlockRel=%e\n", errB_abs, errB_rel);
    printf("errShortLongAbs=%e; errShortLongRel=%e\n", errL_abs, errL_rel);

    free(a);
    free(b);
    free(af);
    free(bf);
    free(a1);
    free(b1);
    free(temp);

    return 0;
}
