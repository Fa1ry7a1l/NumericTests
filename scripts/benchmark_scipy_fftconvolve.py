import argparse
import time

import numpy as np
import scipy
from scipy import signal


def make_input(size: int, rng: np.random.Generator) -> np.ndarray:
    """Generate the same input domain as the C++ benchmark: float values in [0, 1)."""
    return rng.random(size, dtype=np.float64)


def benchmark(signal_size: int, kernel_size: int, runs: int, seed: int, warmup: int) -> None:
    rng = np.random.default_rng(seed)
    x = make_input(signal_size, rng)
    h = make_input(kernel_size, rng)

    result_size = signal_size + kernel_size - 1
    print("SciPy fftconvolve benchmark")
    print(f"NumPy={np.__version__}, SciPy={scipy.__version__}")
    print(
        f"SIGNAL_SIZE={signal_size}, KERNEL_SIZE={kernel_size}, "
        f"RESULT_SIZE={result_size}, RUNS={runs}, WARMUP={warmup}"
    )
    print(f"x.dtype={x.dtype}, h.dtype={h.dtype}")
    print()

    for _ in range(warmup):
        y = signal.fftconvolve(x, h, mode="full")
        if y.size != result_size:
            raise RuntimeError("Unexpected convolution length during warmup")

    times = []
    last = None
    for run in range(runs):
        start = time.perf_counter()
        y = signal.fftconvolve(x, h, mode="full")
        elapsed = time.perf_counter() - start
        if y.size != result_size:
            raise RuntimeError("Unexpected convolution length")
        times.append(elapsed)
        last = y
        print(f"run {run + 1:2d}: scipy.signal.fftconvolve={elapsed:.5f}s")

    arr = np.array(times, dtype=np.float64)
    print()
    print(f"average: {arr.mean():.5f} s")
    print(f"median:  {np.median(arr):.5f} s")
    print(f"min:     {arr.min():.5f} s")
    print(f"max:     {arr.max():.5f} s")
    print(f"std:     {arr.std(ddof=0):.5f} s")
    print(f"last checksum: {float(np.sum(last)):.17g}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signal-size", type=int, default=1_048_576)
    parser.add_argument("--kernel-size", type=int, default=1_024)
    parser.add_argument("--runs", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--warmup", type=int, default=2)
    args = parser.parse_args()

    benchmark(
        signal_size=args.signal_size,
        kernel_size=args.kernel_size,
        runs=args.runs,
        seed=args.seed,
        warmup=args.warmup,
    )


if __name__ == "__main__":
    main()
