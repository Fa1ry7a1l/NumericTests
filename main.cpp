#include <NTL/FFT.h>
#include <NTL/lzz_p.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "third_party/pocketfft/pocketfft_hdronly.h"

#ifndef FIRST_TRANSFORM_LOG2
#define FIRST_TRANSFORM_LOG2 20
#endif

#ifndef LAST_TRANSFORM_LOG2
#define LAST_TRANSFORM_LOG2 24
#endif

#ifndef RUNS
#define RUNS 20
#endif

#ifndef WARMUP
#define WARMUP 3
#endif

#ifndef RANDOM_SEED
#define RANDOM_SEED 42u
#endif

#ifndef NTT_SCALE
#define NTT_SCALE 1000000L
#endif

#ifndef NTT_MODULUS
#define NTT_MODULUS 469762049L
#endif

using Clock = std::chrono::steady_clock;

struct DiffStats {
    long double sse = 0.0L;
    double mse = 0.0;
    double rmse = 0.0;
    double relative_rmse = 0.0;
    double max_abs = 0.0;
};

struct BenchmarkResult {
    size_t size = 0;
    double fft_seconds = 0.0;
    double ntt_core_seconds = 0.0;
    double ntt_full_seconds = 0.0;
    DiffStats fft_diff;
    DiffStats ntt_diff;
};

static constexpr size_t transform_size(long log2_size) {
    return size_t{1} << log2_size;
}

static std::vector<double> make_input(size_t n) {
    std::mt19937_64 rng(RANDOM_SEED);
    std::vector<double> values(n);
    for (double& x : values) {
        x = std::generate_canonical<double, 53>(rng);
    }
    return values;
}

static DiffStats compare_quadratic(const std::vector<double>& reference,
                                   const std::vector<double>& tested) {
    if (reference.size() != tested.size()) {
        throw std::runtime_error("result sizes differ");
    }

    DiffStats stats;
    long double reference_square_sum = 0.0L;
    for (size_t i = 0; i < reference.size(); ++i) {
        const double diff = tested[i] - reference[i];
        stats.sse += static_cast<long double>(diff) * diff;
        reference_square_sum += static_cast<long double>(reference[i]) * reference[i];
        stats.max_abs = std::max(stats.max_abs, std::abs(diff));
    }

    stats.mse = static_cast<double>(stats.sse / static_cast<long double>(reference.size()));
    stats.rmse = std::sqrt(stats.mse);
    stats.relative_rmse = reference_square_sum > 0.0L
                              ? std::sqrt(static_cast<double>(stats.sse / reference_square_sum))
                              : 0.0;
    return stats;
}

static double average(const std::vector<double>& times) {
    return std::accumulate(times.begin(), times.end(), 0.0) / static_cast<double>(times.size());
}

static void pocketfft_forward(const std::vector<double>& input,
                              std::vector<std::complex<double>>& spectrum) {
    const size_t n = input.size();
    pocketfft::shape_t shape{n};
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

static void pocketfft_inverse(const std::vector<std::complex<double>>& spectrum,
                              std::vector<double>& output) {
    const size_t n = output.size();
    pocketfft::shape_t shape{n};
    pocketfft::stride_t real_stride{static_cast<std::ptrdiff_t>(sizeof(double))};
    pocketfft::stride_t complex_stride{static_cast<std::ptrdiff_t>(sizeof(std::complex<double>))};

    pocketfft::c2r(shape,
                   complex_stride,
                   real_stride,
                   0,
                   pocketfft::BACKWARD,
                   spectrum.data(),
                   output.data(),
                   1.0 / static_cast<double>(n),
                   1);
}

static void run_pocketfft_roundtrip(const std::vector<double>& input,
                                    std::vector<double>& work,
                                    std::vector<std::complex<double>>& spectrum,
                                    std::vector<double>& restored) {
    work = input;
    pocketfft_forward(work, spectrum);
    pocketfft_inverse(spectrum, restored);
}

static void scale_to_ntt(const std::vector<double>& input, std::vector<long>& scaled) {
    const long modulus = NTT_MODULUS;
    for (size_t i = 0; i < input.size(); ++i) {
        const long value = static_cast<long>(std::llround(input[i] * static_cast<double>(NTT_SCALE)));
        if (value < 0 || value >= modulus) {
            throw std::runtime_error("NTT_SCALE is too large for NTT_MODULUS");
        }
        scaled[i] = value;
    }
}

static void unscale_from_ntt(const std::vector<long>& recovered, std::vector<double>& restored) {
    const double inv_scale = 1.0 / static_cast<double>(NTT_SCALE);
    for (size_t i = 0; i < recovered.size(); ++i) {
        restored[i] = static_cast<double>(recovered[i]) * inv_scale;
    }
}

static void run_ntl_ntt_core(const NTL::FFTPrimeInfo& prime,
                             long log2_size,
                             const std::vector<long>& scaled,
                             std::vector<long>& spectrum,
                             std::vector<long>& recovered) {
    NTL::FFTFwd(spectrum.data(), scaled.data(), log2_size, prime);
    NTL::FFTRev1(recovered.data(), spectrum.data(), log2_size, prime);
}

static void run_ntl_ntt_full(const NTL::FFTPrimeInfo& prime,
                             long log2_size,
                             const std::vector<double>& input,
                             std::vector<long>& scaled,
                             std::vector<long>& spectrum,
                             std::vector<long>& recovered,
                             std::vector<double>& restored) {
    scale_to_ntt(input, scaled);
    run_ntl_ntt_core(prime, log2_size, scaled, spectrum, recovered);
    unscale_from_ntt(recovered, restored);
}

static BenchmarkResult benchmark_size(long log2_size, const NTL::FFTPrimeInfo& ntt_prime) {
    const size_t n = transform_size(log2_size);
    const size_t spectrum_size = n / 2 + 1;
    const std::vector<double> input = make_input(n);
    std::vector<double> fft_work(n);
    std::vector<double> fft_restored(n);
    std::vector<std::complex<double>> fft_spectrum(spectrum_size);

    std::vector<long> ntt_scaled(n);
    std::vector<long> ntt_spectrum(n);
    std::vector<long> ntt_recovered(n);
    std::vector<double> ntt_restored(n);

    scale_to_ntt(input, ntt_scaled);

    double guard = 0.0;
    for (int i = 0; i < WARMUP; ++i) {
        run_pocketfft_roundtrip(input, fft_work, fft_spectrum, fft_restored);
        run_ntl_ntt_full(ntt_prime,
                         log2_size,
                         input,
                         ntt_scaled,
                         ntt_spectrum,
                         ntt_recovered,
                         ntt_restored);
        guard += fft_restored[static_cast<size_t>(i) % n] + ntt_restored[static_cast<size_t>(i) % n];
    }

    std::vector<double> fft_times;
    std::vector<double> ntt_full_times;
    std::vector<double> ntt_core_times;
    fft_times.reserve(RUNS);
    ntt_full_times.reserve(RUNS);
    ntt_core_times.reserve(RUNS);

    for (int run = 0; run < RUNS; ++run) {
        const auto fft_start = Clock::now();
        run_pocketfft_roundtrip(input, fft_work, fft_spectrum, fft_restored);
        const double fft_elapsed = std::chrono::duration<double>(Clock::now() - fft_start).count();
        fft_times.push_back(fft_elapsed);
        guard += fft_restored[static_cast<size_t>(run * 8191) % n];

        const auto ntt_full_start = Clock::now();
        run_ntl_ntt_full(ntt_prime,
                         log2_size,
                         input,
                         ntt_scaled,
                         ntt_spectrum,
                         ntt_recovered,
                         ntt_restored);
        const double ntt_full_elapsed =
            std::chrono::duration<double>(Clock::now() - ntt_full_start).count();
        ntt_full_times.push_back(ntt_full_elapsed);
        guard += ntt_restored[static_cast<size_t>(run * 6151) % n];

        const auto ntt_core_start = Clock::now();
        run_ntl_ntt_core(ntt_prime, log2_size, ntt_scaled, ntt_spectrum, ntt_recovered);
        const double ntt_core_elapsed =
            std::chrono::duration<double>(Clock::now() - ntt_core_start).count();
        ntt_core_times.push_back(ntt_core_elapsed);
        guard += static_cast<double>(ntt_recovered[static_cast<size_t>(run * 4099) % n]) /
                 static_cast<double>(NTT_SCALE);
    }

    if (guard == -1.0) {
        std::cout << "guard: " << guard << "\n";
    }

    BenchmarkResult result;
    result.size = n;
    result.fft_seconds = average(fft_times);
    result.ntt_core_seconds = average(ntt_core_times);
    result.ntt_full_seconds = average(ntt_full_times);
    result.fft_diff = compare_quadratic(input, fft_restored);
    result.ntt_diff = compare_quadratic(input, ntt_restored);
    return result;
}

static void print_table(const std::vector<BenchmarkResult>& results) {
    std::cout << std::defaultfloat << std::setprecision(6);
    std::cout << "| size | FFT avg, s | NTT core avg, s | speedup FFT/NTT core | NTT with conversion avg, s | speedup FFT/NTT with conversion | NTT RMSE | relative RMSE, % |\n";
    std::cout << "|---:|---:|---:|---:|---:|---:|---:|---:|\n";

    for (const BenchmarkResult& result : results) {
        std::cout << "| " << result.size
                  << " | " << result.fft_seconds
                  << " | " << result.ntt_core_seconds
                  << " | " << (result.fft_seconds / result.ntt_core_seconds)
                  << " | " << result.ntt_full_seconds
                  << " | " << (result.fft_seconds / result.ntt_full_seconds)
                  << " | " << std::scientific << std::setprecision(6) << result.ntt_diff.rmse
                  << " | " << std::defaultfloat << std::setprecision(6)
                  << (result.ntt_diff.relative_rmse * 100.0)
                  << " |\n";
    }
}

int main() {
    static_assert(FIRST_TRANSFORM_LOG2 > 0, "FIRST_TRANSFORM_LOG2 must be positive");
    static_assert(LAST_TRANSFORM_LOG2 >= FIRST_TRANSFORM_LOG2,
                  "LAST_TRANSFORM_LOG2 must be at least FIRST_TRANSFORM_LOG2");
    static_assert(LAST_TRANSFORM_LOG2 <= 26, "469762049 supports NTT lengths up to 2^26");

    if (NTL::CalcMaxRoot(NTT_MODULUS) < LAST_TRANSFORM_LOG2) {
        throw std::runtime_error("NTT_MODULUS does not support requested transform length");
    }

    NTL::zz_p::UserFFTInit(NTT_MODULUS);
    if (NTL::zz_pInfo == nullptr || NTL::zz_pInfo->p_info == nullptr) {
        throw std::runtime_error("NTL did not accept NTT_MODULUS as an FFT prime");
    }
    const NTL::FFTPrimeInfo& ntt_prime = *NTL::zz_pInfo->p_info;

    std::vector<BenchmarkResult> results;
    results.reserve(static_cast<size_t>(LAST_TRANSFORM_LOG2 - FIRST_TRANSFORM_LOG2 + 1));

    std::cout << "FFT/NTT library roundtrip benchmark table\n";
    std::cout << "FIRST_TRANSFORM_LOG2=" << FIRST_TRANSFORM_LOG2
              << ", LAST_TRANSFORM_LOG2=" << LAST_TRANSFORM_LOG2
              << ", RUNS=" << RUNS
              << ", WARMUP=" << WARMUP
              << ", input=[0;1)"
              << ", NTT_MODULUS=" << NTT_MODULUS
              << ", NTT_SCALE=" << NTT_SCALE << "\n\n";

    for (long log2_size = FIRST_TRANSFORM_LOG2; log2_size <= LAST_TRANSFORM_LOG2; ++log2_size) {
        results.push_back(benchmark_size(log2_size, ntt_prime));
    }

    print_table(results);

    return 0;
}
