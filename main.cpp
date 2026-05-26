#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

// Параметры последнего эксперимента.
// SIGNAL_SIZE - длина входного сигнала; KERNEL_SIZE - длина короткого ядра.
#ifndef SIGNAL_SIZE
#define SIGNAL_SIZE 1048576
#endif

#ifndef KERNEL_SIZE
#define KERNEL_SIZE 1024
#endif

#ifndef RUNS
#define RUNS 20
#endif

#ifndef RANDOM_SEED
#define RANDOM_SEED 42u
#endif

// Входные float-данные генерируются в диапазоне [0; 1).
#ifndef INPUT_MIN
#define INPUT_MIN 0.0f
#endif

#ifndef INPUT_MAX
#define INPUT_MAX 1.0f
#endif

// Квантование float -> int16_t: round(x * INTEGER_SCALE).
// Для последнего теста scale=5 дал хороший баланс скорости и точности.
#ifndef INTEGER_SCALE
#define INTEGER_SCALE 5
#endif

// Ускоренный метод: overlap-add по блокам сигнала.
// Каждый блок сворачивается с ядром коротким NTT и добавляется в int64.
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 1024
#endif

// Модуль Ферма 65537 помогает держать короткие NTT; состояние хранится в
// uint32_t, потому что значение 65536 не помещается в uint16_t.
#ifndef SHORT_NTT_MOD
#define SHORT_NTT_MOD 65537u
#endif

#ifndef SHORT_NTT_ROOT
#define SHORT_NTT_ROOT 3u
#endif

using Clock = std::chrono::steady_clock;
using Complex = std::complex<double>;

struct ErrorStats {
    double max_abs = 0.0;
    double mean_abs = 0.0;
    double rmse = 0.0;
    double relative_rmse = 0.0;
    double max_rel = 0.0;
    double mean_rel = 0.0;
};

static size_t next_power_of_two(size_t n) {
    size_t p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

// Эталонный метод: итеративное быстрое преобразование Фурье Cooley-Tukey
// над complex<double>, затем обратное FFT для получения свертки.
static void fft(std::vector<Complex>& a, bool invert) {
    const size_t n = a.size();

    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        const double angle = 2.0 * std::acos(-1.0) / static_cast<double>(len) *
                             (invert ? -1.0 : 1.0);
        const Complex wlen(std::cos(angle), std::sin(angle));

        for (size_t i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for (size_t j = 0; j < len / 2; ++j) {
                const Complex u = a[i + j];
                const Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (Complex& x : a) {
            x /= static_cast<double>(n);
        }
    }
}

static std::vector<double> convolution_fft(const std::vector<float>& signal,
                                           const std::vector<float>& kernel) {
    const size_t result_size = signal.size() + kernel.size() - 1;
    const size_t ntt_size = next_power_of_two(result_size);

    std::vector<Complex> a(ntt_size);
    std::vector<Complex> b(ntt_size);

    for (size_t i = 0; i < signal.size(); ++i) {
        a[i] = Complex(signal[i], 0.0);
    }
    for (size_t i = 0; i < kernel.size(); ++i) {
        b[i] = Complex(kernel[i], 0.0);
    }

    fft(a, false);
    fft(b, false);
    for (size_t i = 0; i < ntt_size; ++i) {
        a[i] *= b[i];
    }
    fft(a, true);

    std::vector<double> result(result_size);
    for (size_t i = 0; i < result_size; ++i) {
        result[i] = a[i].real();
    }
    return result;
}

static uint32_t mod_pow(uint32_t base, uint32_t exp, uint32_t mod) {
    uint64_t result = 1;
    uint64_t x = base;
    while (exp > 0) {
        if (exp & 1u) {
            result = (result * x) % mod;
        }
        x = (x * x) % mod;
        exp >>= 1u;
    }
    return static_cast<uint32_t>(result);
}

// Быстрое приведение по модулю 65537: так как 65536 == -1 (mod 65537),
// x mod 65537 можно получить через разность младших и старших 16 бит.
static uint32_t reduce_fermat_65537(uint64_t x) {
    int64_t reduced = static_cast<int64_t>(x & 0xffffu) - static_cast<int64_t>(x >> 16u);
    if (reduced < 0) {
        reduced += SHORT_NTT_MOD;
    }
    if (reduced >= static_cast<int64_t>(SHORT_NTT_MOD)) {
        reduced -= SHORT_NTT_MOD;
    }
    return static_cast<uint32_t>(reduced);
}

static uint32_t mod_mul(uint32_t a, uint32_t b) {
    return reduce_fermat_65537(static_cast<uint64_t>(a) * static_cast<uint64_t>(b));
}

// План короткого NTT: bit-reversal и корни единицы считаются один раз.
// Это аналог FFT-plan: сам алгоритм тот же Cooley-Tukey, но без лишней
// подготовки на каждом блоке.
struct ShortNttPlan {
    size_t n = 0;
    std::vector<size_t> bit_reverse;
    std::vector<uint32_t> forward_wlen;
    std::vector<uint32_t> inverse_wlen;
    uint32_t inverse_n = 1;

    explicit ShortNttPlan(size_t transform_size) : n(transform_size), bit_reverse(transform_size) {
        if (n == 0 || (n & (n - 1)) != 0) {
            throw std::runtime_error("NTT size must be a power of two");
        }
        if (n > SHORT_NTT_MOD - 1u) {
            throw std::runtime_error("Short NTT does not support this block size");
        }

        for (size_t i = 1, j = 0; i < n; ++i) {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;
            bit_reverse[i] = j;
        }

        for (size_t len = 2; len <= n; len <<= 1) {
            const uint32_t wlen =
                mod_pow(SHORT_NTT_ROOT, (SHORT_NTT_MOD - 1u) / static_cast<uint32_t>(len), SHORT_NTT_MOD);
            const uint32_t inv_wlen = mod_pow(wlen, SHORT_NTT_MOD - 2u, SHORT_NTT_MOD);
            forward_wlen.push_back(wlen);
            inverse_wlen.push_back(inv_wlen);
        }

        inverse_n = mod_pow(static_cast<uint32_t>(n), SHORT_NTT_MOD - 2u, SHORT_NTT_MOD);
    }

    void apply(std::vector<uint32_t>& a, bool invert) const {
        for (size_t i = 1; i < n; ++i) {
            const size_t j = bit_reverse[i];
            if (i < j) {
                std::swap(a[i], a[j]);
            }
        }

        const std::vector<uint32_t>& wlen_by_stage = invert ? inverse_wlen : forward_wlen;
        size_t stage = 0;
        for (size_t len = 2; len <= n; len <<= 1, ++stage) {
            const size_t half = len / 2;
            const uint32_t wlen = wlen_by_stage[stage];
            for (size_t i = 0; i < n; i += len) {
                uint32_t w = 1;
                for (size_t j = 0; j < half; ++j) {
                    const uint32_t u = a[i + j];
                    const uint32_t v = mod_mul(a[i + j + half], w);
                    const uint32_t sum = u + v;
                    a[i + j] = (sum >= SHORT_NTT_MOD) ? (sum - SHORT_NTT_MOD) : sum;
                    a[i + j + half] = (u >= v) ? (u - v) : (u + SHORT_NTT_MOD - v);
                    w = mod_mul(w, wlen);
                }
            }
        }

        if (invert) {
            for (uint32_t& x : a) {
                x = mod_mul(x, inverse_n);
            }
        }
    }
};

static int16_t quantize_to_int16(float x) {
    const double scaled = static_cast<double>(x) * static_cast<double>(INTEGER_SCALE);
    const long long rounded = std::llround(scaled);
    if (rounded < std::numeric_limits<int16_t>::min() ||
        rounded > std::numeric_limits<int16_t>::max()) {
        throw std::runtime_error("INTEGER_SCALE does not fit into int16_t for this input");
    }
    return static_cast<int16_t>(rounded);
}

static uint32_t encode_mod(int16_t value) {
    const int signed_value = static_cast<int>(value);
    if (signed_value >= 0) {
        return static_cast<uint32_t>(signed_value);
    }
    return SHORT_NTT_MOD - static_cast<uint32_t>(-signed_value);
}

static int64_t decode_mod(uint32_t value) {
    if (value > SHORT_NTT_MOD / 2u) {
        return static_cast<int64_t>(value) - static_cast<int64_t>(SHORT_NTT_MOD);
    }
    return static_cast<int64_t>(value);
}

static void add_short_ntt_block(const int16_t* block,
                                size_t block_len,
                                const std::vector<uint32_t>& kernel_ntt,
                                size_t kernel_size,
                                const ShortNttPlan& ntt_plan,
                                std::vector<uint32_t>& work,
                                std::vector<int32_t>& accumulated,
                                size_t offset) {
    const size_t result_size = block_len + kernel_size - 1;
    const size_t transform_size = work.size();

    std::fill(work.begin(), work.end(), 0u);
    for (size_t i = 0; i < block_len; ++i) {
        work[i] = encode_mod(block[i]);
    }

    ntt_plan.apply(work, false);
    for (size_t i = 0; i < transform_size; ++i) {
        work[i] = mod_mul(work[i], kernel_ntt[i]);
    }
    ntt_plan.apply(work, true);

    for (size_t i = 0; i < result_size; ++i) {
        accumulated[offset + i] += static_cast<int32_t>(decode_mod(work[i]));
    }
}

struct BlockedNttPlan {
    size_t signal_size = 0;
    size_t kernel_size = 0;
    size_t result_size = 0;
    size_t transform_size = 0;
    std::vector<int16_t> signal_i16;
    std::vector<uint32_t> kernel_ntt;
    ShortNttPlan ntt_plan;

    BlockedNttPlan(const std::vector<float>& signal, const std::vector<float>& kernel)
        : signal_size(signal.size()),
          kernel_size(kernel.size()),
          result_size(signal.size() + kernel.size() - 1),
          transform_size(next_power_of_two(BLOCK_SIZE + kernel.size() - 1)),
          signal_i16(signal.size()),
          kernel_ntt(transform_size),
          ntt_plan(transform_size) {
        std::vector<int16_t> kernel_i16(kernel.size());
        std::transform(signal.begin(), signal.end(), signal_i16.begin(), quantize_to_int16);
        std::transform(kernel.begin(), kernel.end(), kernel_i16.begin(), quantize_to_int16);

        // Kernel один и тот же для всех запусков и всех блоков, поэтому его
        // NTT-спектр входит в подготовленный план.
        for (size_t i = 0; i < kernel_i16.size(); ++i) {
            kernel_ntt[i] = encode_mod(kernel_i16[i]);
        }
        ntt_plan.apply(kernel_ntt, false);
    }
};

struct BlockedNttWorkspace {
    std::vector<int32_t> accumulated;
    std::vector<uint32_t> work;
    std::vector<double> result;

    explicit BlockedNttWorkspace(const BlockedNttPlan& plan)
        : accumulated(plan.result_size), work(plan.transform_size), result(plan.result_size) {}
};

// Ускоренный код сравнения: подготовленный float -> int16_t план выполняет
// overlap-add из коротких NTT. Частичные результаты блоков накапливаются
// в int32_t, затем возвращаются в double делением на INTEGER_SCALE^2.
static const std::vector<double>& execute_blocked_int16_ntt(const BlockedNttPlan& plan,
                                                            BlockedNttWorkspace& workspace,
                                                            double& elapsed_seconds) {
    const auto start = Clock::now();

    std::fill(workspace.accumulated.begin(), workspace.accumulated.end(), 0);

    for (size_t offset = 0; offset < plan.signal_i16.size(); offset += BLOCK_SIZE) {
        const size_t block_len = std::min<size_t>(BLOCK_SIZE, plan.signal_i16.size() - offset);
        add_short_ntt_block(plan.signal_i16.data() + offset,
                            block_len,
                            plan.kernel_ntt,
                            plan.kernel_size,
                            plan.ntt_plan,
                            workspace.work,
                            workspace.accumulated,
                            offset);
    }

    const double scale_back = static_cast<double>(INTEGER_SCALE) * static_cast<double>(INTEGER_SCALE);
    for (size_t i = 0; i < plan.result_size; ++i) {
        workspace.result[i] = static_cast<double>(workspace.accumulated[i]) / scale_back;
    }

    elapsed_seconds = std::chrono::duration<double>(Clock::now() - start).count();
    return workspace.result;
}

static ErrorStats compare_results(const std::vector<double>& reference,
                                  const std::vector<double>& tested) {
    if (reference.size() != tested.size()) {
        throw std::runtime_error("Result sizes differ");
    }

    ErrorStats stats;
    long double abs_sum = 0.0;
    long double sq_sum = 0.0;
    long double ref_sq_sum = 0.0;
    long double rel_sum = 0.0;
    size_t rel_count = 0;

    for (size_t i = 0; i < reference.size(); ++i) {
        const double diff = tested[i] - reference[i];
        const double abs_diff = std::abs(diff);
        stats.max_abs = std::max(stats.max_abs, abs_diff);
        abs_sum += abs_diff;
        sq_sum += static_cast<long double>(diff) * diff;
        ref_sq_sum += static_cast<long double>(reference[i]) * reference[i];

        const double denom = std::abs(reference[i]);
        if (denom > 1e-12) {
            const double rel = abs_diff / denom;
            stats.max_rel = std::max(stats.max_rel, rel);
            rel_sum += rel;
            ++rel_count;
        }
    }

    const long double n = static_cast<long double>(reference.size());
    stats.mean_abs = static_cast<double>(abs_sum / n);
    stats.rmse = std::sqrt(static_cast<double>(sq_sum / n));
    stats.relative_rmse = (ref_sq_sum > 0.0L)
                              ? std::sqrt(static_cast<double>(sq_sum / ref_sq_sum))
                              : 0.0;
    stats.mean_rel = rel_count ? static_cast<double>(rel_sum / rel_count) : 0.0;
    return stats;
}

static std::vector<float> make_input(size_t n, std::mt19937& rng) {
    std::vector<float> values(n);
    for (float& x : values) {
        const float unit = std::generate_canonical<float, 24>(rng);
        x = INPUT_MIN + (INPUT_MAX - INPUT_MIN) * unit;
    }
    return values;
}

int main() {
    static_assert(SHORT_NTT_MOD == 65537u,
                  "Fast modular reduction is specialized for modulus 65537");
    static_assert(INPUT_MIN == 0.0f && INPUT_MAX == 1.0f,
                  "The last test expects float input in [0; 1)");
    static_assert(BLOCK_SIZE + KERNEL_SIZE - 1 <= SHORT_NTT_MOD - 1,
                  "For mod 65537, the short NTT size must not exceed 65536");

    std::mt19937 rng(RANDOM_SEED);
    const std::vector<float> signal = make_input(SIGNAL_SIZE, rng);
    const std::vector<float> kernel = make_input(KERNEL_SIZE, rng);

    const auto prepare_start = Clock::now();
    const BlockedNttPlan accelerated_plan(signal, kernel);
    const double prepare_seconds = std::chrono::duration<double>(Clock::now() - prepare_start).count();
    BlockedNttWorkspace accelerated_workspace(accelerated_plan);

    double fft_total = 0.0;
    double accelerated_total = 0.0;
    ErrorStats last_stats;

    std::cout << "Last benchmark only\n";
    std::cout << "SIGNAL_SIZE=" << SIGNAL_SIZE
              << ", KERNEL_SIZE=" << KERNEL_SIZE
              << ", RUNS=" << RUNS
              << ", INTEGER_SCALE=" << INTEGER_SCALE
              << ", BLOCK_SIZE=" << BLOCK_SIZE
              << ", NTT_MOD=" << SHORT_NTT_MOD << "\n\n";

    for (int run = 0; run < RUNS; ++run) {
        const auto fft_start = Clock::now();
        const std::vector<double> reference = convolution_fft(signal, kernel);
        const auto fft_end = Clock::now();
        fft_total += std::chrono::duration<double>(fft_end - fft_start).count();

        double accelerated_seconds = 0.0;
        const std::vector<double>& accelerated =
            execute_blocked_int16_ntt(accelerated_plan, accelerated_workspace, accelerated_seconds);

        accelerated_total += accelerated_seconds;
        last_stats = compare_results(reference, accelerated);

        std::cout << "run " << std::setw(2) << (run + 1)
                  << ": FFT=" << std::fixed << std::setprecision(5)
                  << std::chrono::duration<double>(fft_end - fft_start).count()
                  << "s, prepared int16+NTT=" << accelerated_seconds
                  << "s\n";
    }

    const double fft_avg = fft_total / RUNS;
    const double accelerated_avg = accelerated_total / RUNS;
    const double accelerated_with_amortized_prepare = accelerated_avg + prepare_seconds / RUNS;

    std::cout << "\nAverage over " << RUNS << " runs\n";
    std::cout << "FFT Cooley-Tukey reference:          " << fft_avg << " s\n";
    std::cout << "one-time int16+NTT plan preparation: " << prepare_seconds << " s\n";
    std::cout << "prepared blocked int16 NTT execute:  " << accelerated_avg << " s\n";
    std::cout << "execute + amortized preparation:     " << accelerated_with_amortized_prepare << " s\n";
    std::cout << "speedup vs FFT, prepared execute:    " << (fft_avg / accelerated_avg) << "x\n";
    std::cout << "speedup vs FFT, amortized prepare:   " << (fft_avg / accelerated_with_amortized_prepare) << "x\n\n";

    std::cout << "Accuracy vs FFT reference\n";
    std::cout << "max abs:        " << last_stats.max_abs << "\n";
    std::cout << "mean abs:       " << last_stats.mean_abs << "\n";
    std::cout << "RMSE:           " << last_stats.rmse << "\n";
    std::cout << "relative RMSE:  " << (last_stats.relative_rmse * 100.0) << "%\n";
    std::cout << "max rel:        " << (last_stats.max_rel * 100.0) << "%\n";
    std::cout << "mean rel:       " << (last_stats.mean_rel * 100.0) << "%\n";

    return 0;
}
