#include <algorithm>
#include <chrono>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "../third_party/pocketfft/pocketfft_hdronly.h"

#ifndef SIGNAL_SIZE
#define SIGNAL_SIZE 1048576
#endif

#ifndef KERNEL_SIZE
#define KERNEL_SIZE 1024
#endif

#ifndef RUNS
#define RUNS 20
#endif

#ifndef WARMUP
#define WARMUP 2
#endif

#ifndef RANDOM_SEED
#define RANDOM_SEED 42u
#endif

using Clock = std::chrono::steady_clock;

static std::vector<double> make_input(size_t n, std::mt19937& rng) {
    std::vector<double> values(n);
    for (double& x : values) {
        x = std::generate_canonical<double, 53>(rng);
    }
    return values;
}

static void print_stats(const std::vector<double>& times) {
    std::vector<double> sorted = times;
    std::sort(sorted.begin(), sorted.end());
    const double sum = std::accumulate(times.begin(), times.end(), 0.0);
    const double mean = sum / static_cast<double>(times.size());
    const double median = sorted[sorted.size() / 2];
    double sq = 0.0;
    for (double t : times) {
        const double d = t - mean;
        sq += d * d;
    }

    std::cout << "average: " << mean << " s\n";
    std::cout << "median:  " << median << " s\n";
    std::cout << "min:     " << sorted.front() << " s\n";
    std::cout << "max:     " << sorted.back() << " s\n";
    std::cout << "std:     " << std::sqrt(sq / static_cast<double>(times.size())) << " s\n";
}

static void real_fft_forward(const std::vector<double>& input,
                             std::vector<std::complex<double>>& spectrum,
                             size_t nfft) {
    pocketfft::shape_t shape{nfft};
    pocketfft::stride_t real_stride{static_cast<std::ptrdiff_t>(sizeof(double))};
    pocketfft::stride_t complex_stride{static_cast<std::ptrdiff_t>(sizeof(std::complex<double>))};
    pocketfft::r2c(shape,
                   real_stride,
                   complex_stride,
                   0,
                   pocketfft::FORWARD,
                   input.data(),
                   spectrum.data(),
                   1.0,
                   1);
}

static void real_fft_inverse(const std::vector<std::complex<double>>& spectrum,
                             std::vector<double>& output,
                             size_t nfft) {
    pocketfft::shape_t shape{nfft};
    pocketfft::stride_t real_stride{static_cast<std::ptrdiff_t>(sizeof(double))};
    pocketfft::stride_t complex_stride{static_cast<std::ptrdiff_t>(sizeof(std::complex<double>))};
    pocketfft::c2r(shape,
                   complex_stride,
                   real_stride,
                   0,
                   pocketfft::BACKWARD,
                   spectrum.data(),
                   output.data(),
                   1.0 / static_cast<double>(nfft),
                   1);
}

static double checksum(const std::vector<double>& data, size_t n) {
    return std::accumulate(data.begin(), data.begin() + static_cast<std::ptrdiff_t>(n), 0.0);
}

int main() {
    const size_t signal_size = SIGNAL_SIZE;
    const size_t kernel_size = KERNEL_SIZE;
    const size_t result_size = signal_size + kernel_size - 1;

    // PocketFFT can use FFTPACK-style mixed-radix sizes. This mirrors the idea
    // of scipy.fft.next_fast_len better than forcing the next power of two.
    const size_t nfft = pocketfft::detail::util::good_size_real(result_size);
    const size_t spectrum_size = nfft / 2 + 1;

    std::mt19937 rng(RANDOM_SEED);
    const std::vector<double> signal = make_input(signal_size, rng);
    const std::vector<double> kernel = make_input(kernel_size, rng);

    std::vector<double> padded_signal(nfft, 0.0);
    std::vector<double> padded_kernel(nfft, 0.0);
    std::copy(signal.begin(), signal.end(), padded_signal.begin());
    std::copy(kernel.begin(), kernel.end(), padded_kernel.begin());

    std::vector<double> x(nfft);
    std::vector<double> h(nfft);
    std::vector<double> y(nfft);
    std::vector<std::complex<double>> x_freq(spectrum_size);
    std::vector<std::complex<double>> h_freq(spectrum_size);
    std::vector<std::complex<double>> product(spectrum_size);

    std::cout << "PocketFFT C++ benchmark\n";
    std::cout << "SIGNAL_SIZE=" << signal_size
              << ", KERNEL_SIZE=" << kernel_size
              << ", RESULT_SIZE=" << result_size
              << ", NFFT=" << nfft
              << ", RUNS=" << RUNS
              << ", WARMUP=" << WARMUP << "\n\n";

    for (int i = 0; i < WARMUP; ++i) {
        x = padded_signal;
        h = padded_kernel;
        real_fft_forward(x, x_freq, nfft);
        real_fft_forward(h, h_freq, nfft);
        for (size_t j = 0; j < spectrum_size; ++j) {
            product[j] = x_freq[j] * h_freq[j];
        }
        real_fft_inverse(product, y, nfft);
    }

    std::vector<double> full_times;
    full_times.reserve(RUNS);
    for (int run = 0; run < RUNS; ++run) {
        const auto start = Clock::now();
        x = padded_signal;
        h = padded_kernel;
        real_fft_forward(x, x_freq, nfft);
        real_fft_forward(h, h_freq, nfft);
        for (size_t j = 0; j < spectrum_size; ++j) {
            product[j] = x_freq[j] * h_freq[j];
        }
        real_fft_inverse(product, y, nfft);
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        full_times.push_back(elapsed);
        std::cout << "full run " << std::setw(2) << (run + 1)
                  << ": pocketfft fftconvolve=" << std::fixed << std::setprecision(5)
                  << elapsed << "s\n";
    }

    std::cout << "\nPocketFFT full convolution, kernel FFT inside every run\n";
    print_stats(full_times);
    std::cout << "last checksum: " << std::setprecision(17) << checksum(y, result_size) << "\n\n";

    h = padded_kernel;
    const auto prepare_start = Clock::now();
    real_fft_forward(h, h_freq, nfft);
    const double prepare_seconds = std::chrono::duration<double>(Clock::now() - prepare_start).count();

    for (int i = 0; i < WARMUP; ++i) {
        x = padded_signal;
        real_fft_forward(x, x_freq, nfft);
        for (size_t j = 0; j < spectrum_size; ++j) {
            product[j] = x_freq[j] * h_freq[j];
        }
        real_fft_inverse(product, y, nfft);
    }

    std::vector<double> prepared_times;
    prepared_times.reserve(RUNS);
    for (int run = 0; run < RUNS; ++run) {
        const auto start = Clock::now();
        x = padded_signal;
        real_fft_forward(x, x_freq, nfft);
        for (size_t j = 0; j < spectrum_size; ++j) {
            product[j] = x_freq[j] * h_freq[j];
        }
        real_fft_inverse(product, y, nfft);
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        prepared_times.push_back(elapsed);
        std::cout << "prepared run " << std::setw(2) << (run + 1)
                  << ": pocketfft prepared-kernel=" << std::fixed << std::setprecision(5)
                  << elapsed << "s\n";
    }

    std::cout << "\nPocketFFT prepared kernel convolution\n";
    std::cout << "one-time kernel FFT preparation: " << prepare_seconds << " s\n";
    print_stats(prepared_times);
    std::cout << "last checksum: " << std::setprecision(17) << checksum(y, result_size) << "\n";

    return 0;
}
