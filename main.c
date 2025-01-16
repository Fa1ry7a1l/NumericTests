#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

int n = 128;
int M = 2048;
int M1 = 2048;
int M2 = 2048;
double** a;
double* b;
double* x;
double* res;
double* res2;
int iterationsCount = 10000;
int epochs = 100;


long long time_spent_usual = 0;
long long time_spent_optimized = 0;

void generateMatrix(double** a)
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < n; j++)
        {
            a[i][j] = rand() / (double)((long)RAND_MAX + 1);
            sum += a[i][j];
        }
        sum += 0.03 * sum;

        for (int j = 0; j < n; j++) //todo не должно быть нормальзации
        {
            a[i][j] /= sum;
        }
    }
}

void free_ram(void)
{
    for (int i = 0; i < n; i++)
    {
        free(a[i]);
    }
    free(b);
    free(res);
    free(x);
}

void calculate_floats(void)
{
    double* x_float = (double*)malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++)
    {
        x_float[i] = x[i];
    }

    time_t begin = clock();

    for (int it = 0; it < iterationsCount; it++)
    {
        for (int i = 0; i < n; i++)
        {
            res[i] = 0;
            for (int j = 0; j < n; j++)
            {
                res[i] += a[i][j] * x_float[j];
            }
            res[i] += b[i];
        }

        for (int i = 0; i < n; i++)
        {
            x_float[i] = res[i];
        }
    }
    time_spent_usual = clock() - begin;

    free(x_float);
}

void calculate_ints(void)
{
    unsigned int* x_int = (unsigned int*)malloc(sizeof(unsigned int) * n);
    for (int i = 0; i < n; i++) //домножаем вектор x на M
    {
        x_int[i] = x[i] * M;
    }
    unsigned int** a_int = (unsigned int**)malloc(sizeof(unsigned int*) * n);
    for (int i = 0; i < n; i++)//домножаем матрицу а на M
    {
        a_int[i] = (unsigned int*)malloc(sizeof(unsigned int) * n);
        for (int j = 0; j < n; j++)
        {
            a_int[i][j] = a[i][j] * M;
        }
    }
    unsigned int* b_int = (unsigned int*)malloc(sizeof(unsigned int) * n);
    for (int i = 0; i < n; i++)//домножаем вектор b на M*М
    {
        b_int[i] = b[i] * M * M;
    }
    unsigned int* res_int = (unsigned int*)malloc(sizeof(unsigned int) * n);
    for (int i = 0; i < n; i++)
    {
        res_int[i] = 0;
    }

    time_t begin = clock();

    for (int it = 0; it < iterationsCount; it++)
    {
        for (int i = 0; i < n; i++)
        {
            int temp = 0;
            for (int j = 0; j < n; j++)
            {
                temp += a_int[i][j] * x_int[j];
            }
            temp += b_int[i]; //посчитали в целых числах перемножение строки матрицы на вектор x + соответствующий элемент b
            res2[i] = temp / M; //делим на M
        }

        for (int i = 0; i < n; i++)
        {
            x_int[i] = res2[i]; //переприсываиваем результаты вычисления во входной вектор, чтобы пойти на новую итерацию.
        }
    }

    time_spent_optimized = clock() - begin;

    for (int i = 0; i < n; i++)
    {
        res2[i] = res2[i] / (double)M;//делим результирующий вектор еще раз на M
    }


    free(x_int);
    for (int i = 0; i < n; i++)
    {
        free(a_int[i]);
    }
    free(a_int);
    free(b_int);
    free(res_int);
}

void print_errors(void)
{
    double delta1 = 0, delta2 = 0, norm1 = 0, norm2 = 0, relativeError1 = 0, relativeError2 = 0;

    for (int i = 0; i < n; i++)
    {
        delta1 += fabs(res[i] - res2[i]);
        norm1 += fabs(res[i]);

        delta2 += pow(res[i] - res2[i], 2);
        norm2 += pow(res[i], 2);
    }

    relativeError1 = delta1 * (double)100 / norm1;
    relativeError2 = sqrt(delta2) * (double)100 / sqrt(norm2);

    printf("Relative Error1: %f norm: %f\n", relativeError1, norm1);
    printf("Relative Error2: %f norm: %f\n", relativeError2, norm2);
}

void print_not_null(void)
{
    int a_count = 0, b_count = 0, x_count = 0, y_count = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (a[i][j] > FLT_EPSILON || a[i][j] < -FLT_EPSILON)
            {
                // printf("A[%d][%d]: %f\n", i, j, a[i][j]);
                a_count++;
            }
        }
        if (b[i] > FLT_EPSILON || b[i] < -FLT_EPSILON)
        {
            b_count++;
        }

        if (x[i] > FLT_EPSILON || x[i] < -FLT_EPSILON)
        {
            x_count++;
        }
        if (res[i] > FLT_EPSILON || res[i] < -FLT_EPSILON)
        {
            y_count++;
        }
    }

    printf("a_count: %d, b_count: %d, x_count: %d, y_count: %d\n", a_count, b_count, x_count, y_count);
}

int main(void)
{
    long long seed = time(NULL);
    printf("seed: %lld\n", seed);
    srand(seed);


    a = (double**)malloc(sizeof(double*) * n);
    for (int i = 0; i < n; i++)
    {
        a[i] = (double*)malloc(sizeof(double) * n);
        for (int j = 0; j < n; j++)
        {
            a[i][j] = 0;
        }
    }
    generateMatrix(a);
    b = (double*)malloc(sizeof(double) * n);
    x = (double*)malloc(sizeof(double) * n);
    res = (double*)malloc(sizeof(double) * n);
    res2 = (double*)malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++)
    {
        b[i] = rand() / (double)((long)RAND_MAX + 1);
        x[i] = rand() / (double)((long)RAND_MAX + 1);
        res[i] = 0;
        res2[i] = 0;
    }

    for (int i = 0; i < epochs; i++)
    {
        calculate_floats();
        calculate_ints();
        if (i % 10 == 0)
            printf("IterationsCount: %d\n", i);
    }

    print_not_null();

    print_errors();

    printf("some noize1%f\n", res[0]);
    printf("some noize2%f\n", res2[0]);
    printf("Time spent optimized: %f\n", (double)time_spent_optimized / CLOCKS_PER_SEC / epochs);
    printf("Time spent usual: %f\n", (double)time_spent_usual / CLOCKS_PER_SEC / epochs);

    free_ram();
    return 0;
}
