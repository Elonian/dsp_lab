#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

double complex *NAIVE_FFT(double A[], int size)
{
    double complex **F = (double complex **)malloc(size * sizeof(double complex *));
    for (int i = 0; i < size; i++)
    {
        F[i] = (double complex *)malloc(size * sizeof(double complex));
    }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            F[i][j] = cexp(((-1) * I * 2 * M_PI * i * j) / size);
        }
    }

    double complex *ans = (double complex *)malloc(size * sizeof(double complex));
    for (int i = 0; i < size; i++)
    {
        double complex temp = 0;
        int j = 0;
        for (j = 0; j < size; j++)
        {
            temp = temp + F[i][j] * A[j];
        }
        ans[i] = temp;
    }
    free(F);
    return ans;
}

double complex *NAIVE_IFFT(double complex A[], int size)
{
    double complex F[size][size];
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            F[i][j] = cexp((I * 2 * M_PI * i * j) / size) / size;
        }
    }
    double complex *ans = (double complex *)malloc(size * sizeof(double complex));
    for (int i = 0; i < size; i++)
    {
        double complex temp = 0;
        int j = 0;
        for (j = 0; j < size; j++)
        {
            temp = temp + F[i][j] * A[j];
        }
        ans[i] = temp;
    }
    return ans;
}

void FFT(double complex *x, int data_size)
{
    if (data_size == 1)
    {
        return;
    }
    double complex *odd = (double complex *)malloc(data_size / 2 * sizeof(double complex));
    double complex *even = (double complex *)malloc(data_size / 2 * sizeof(double complex));
    for (int i = 0; i < data_size / 2; i++)
    {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    FFT(even, data_size / 2);
    FFT(odd, data_size / 2);

    for (int i = 0; i < data_size / 2; i++)
    {
        double complex e = CMPLX(cos(2 * M_PI * i / data_size), -sin(2 * M_PI * i / data_size));

        x[i] = even[i] + e * odd[i];
        x[i + data_size / 2] = even[i] - e * odd[i];
    }
    free(even);
    free(odd);
}

double complex polyval(double coeff[], int size, double complex x)
{
    double complex sum = coeff[0];
    for (int i = 1; i < size; i++)
    {
        sum = sum * (x) + coeff[i];
    }
    return sum;
}

double complex *H(double num[], double den[], int num_size, double complex z[], int size_data)
{
    double complex *n = (double complex *)malloc(size_data * sizeof(double complex));
    for (int i = 0; i < size_data; i++)
    {
        n[i] = polyval(num, num_size, z[i]);
    }
    double complex *d = (double complex *)malloc(size_data * sizeof(double complex));
    for (int i = 0; i < size_data; i++)
    {
        d[i] = polyval(den, num_size, z[i]);
    }

    double complex *ans = (double complex *)malloc(size_data * sizeof(double complex));

    for (int i = 0; i < size_data; i++)
    {
        ans[i] = n[i] / d[i];
    }
    return ans;
}

void IFFT(double complex *x, int data_size)
{
    if (data_size == 1)
    {
        return;
    }

    double complex *even = (double complex *)malloc(data_size / 2 * sizeof(double complex));
    double complex *odd = (double complex *)malloc(data_size / 2 * sizeof(double complex));
    for (int i = 0; i < data_size / 2; i++)
    {
        even[i] = x[2 * i];
        odd[i] = x[(2 * i) + 1];
    }

    IFFT(even, data_size / 2);
    IFFT(odd, data_size / 2);

    for (int i = 0; 2 * i < data_size; i++)
    {
        double complex e = CMPLX(cos(2 * M_PI * i / data_size), sin(2 * M_PI * i / data_size));

        x[i] = even[i] + e * odd[i];
        x[i + data_size / 2] = even[i] - e * odd[i];
    }
    free(even);
    free(odd);
}

void my_filter(double complex *x, double complex *H, long int n)
{
    FILE *fouty, *foutYFFT;
    FFT(x, n);
    double complex *Y = malloc(n * sizeof(double complex));
    for (int i = 0; i < n; i++)
    {
        Y[i] = H[i] * x[i];
    }
    foutYFFT = fopen("../data/YFFT.dat", "w");
    for (int i = 0; i < n; i++)
    {
        fprintf(foutYFFT, "%lf+j%lf\n", creal(Y[i]), cimag(Y[i]));
    }
    fclose(foutYFFT);
    IFFT(Y, n);
    fouty = fopen("../data/y.dat", "w");
    for (int i = 0; i < n; i++)
    {
        Y[i] = Y[i] / n;
        fprintf(fouty, "%lf %lf\n", creal(Y[i]), cimag(Y[i]));
    }
    fclose(fouty);
}

int main()
{
    long int n = pow(2, 21);
    double complex *x = (double complex *)malloc(n * sizeof(double complex));
    FILE *finx, *finH;
    finx = fopen("../data/x.dat", "r");

    int k = 0;
    while (!feof(finx) && k < n)
    {
        double val;
        fscanf(finx, "%lf", &val);
        k++;
        x[k] = CMPLX(val, 0);
    }
    fclose(finx);

    double complex *H = malloc(n * sizeof(double complex));
    finH = fopen("../data/H_z.dat", "r");
    double r, i;
    k = 0;
    while (!feof(finH) && k < n)
    {
        fscanf(finH, "%lf %lf", &r, &i);
        k++;
        H[k] = CMPLX(r, i);
    }
    my_filter(x, H, n);

    // printf("ok");
    return 0;
}
