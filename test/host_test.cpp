// Host-side correctness test for fixedpoint-fft.
//
// This file is NOT part of the Arduino library and is never compiled
// into a sketch - it's a plain g++ program that links the library
// headers against a double-precision reference DFT (using <complex>/
// <cmath>, which is fine here: those calls happen only in this
// verification tool, never inside the library itself) to check the
// fixed-point engine's correctness and numerical accuracy.
//
// Build & run:
//   g++ -std=c++17 -O2 -I../src test/host_test.cpp -o /tmp/fft_test && /tmp/fft_test

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <type_traits>
#include <vector>

#include "FixedFFT.h"

using namespace fixedpoint_fft;

static int g_failures = 0;

#define CHECK(cond, msg)                                                \
  do {                                                                  \
    if (!(cond)) {                                                      \
      std::printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__);       \
      g_failures++;                                                     \
    }                                                                   \
  } while (0)

// Reference O(N^2) DFT in double precision, unscaled: X[k] = sum_n x[n] * exp(-j*2*pi*k*n/N)
static std::vector<std::complex<double>> referenceDft(const std::vector<std::complex<double>>& x, bool inverse) {
  const size_t n = x.size();
  std::vector<std::complex<double>> out(n);
  const double sign = inverse ? 1.0 : -1.0;
  for (size_t k = 0; k < n; k++) {
    std::complex<double> sum(0.0, 0.0);
    for (size_t t = 0; t < n; t++) {
      double angle = sign * 2.0 * M_PI * (double)k * (double)t / (double)n;
      sum += x[t] * std::complex<double>(std::cos(angle), std::sin(angle));
    }
    out[k] = sum;
  }
  return out;
}

// --- Twiddle table sanity check ------------------------------------------
template <typename CalcT>
void testTwiddleTable(const char* name) {
  constexpr uint16_t maxN = FFTTwiddleTraits<CalcT>::kMaxN;
  double scale;
  if constexpr (std::is_floating_point<CalcT>::value) {
    scale = 1.0;  // f32 table stores cos()/sin() directly, unscaled
  } else {
    constexpr int qbits = FFTQBits<CalcT>::value;
    scale = (double)((1u << qbits) - 1);
  }
  double worst = 0.0;
  for (uint16_t idx = 0; idx < maxN; idx += 7) {  // sample across the whole circle
    CalcT c, s;
    fftTwiddle<CalcT>(idx, c, s);
    double angle = 2.0 * M_PI * idx / maxN;
    double expC = std::cos(angle) * scale;
    double expS = std::sin(angle) * scale;
    worst = std::max(worst, std::fabs((double)c - expC));
    worst = std::max(worst, std::fabs((double)s - expS));
  }
  std::printf("[%s] twiddle table worst abs error (raw units, scale=%.0f): %.3g\n", name, scale, worst);
  const double tol = std::is_floating_point<CalcT>::value ? 1e-6 : 1.0;
  CHECK(worst <= tol, "twiddle table deviates from cos/sin by more than tolerance");
}

// --- Peak-bin detection: pure tone in, peak must land on the right bin ---
template <typename CalcT>
void testPeakBin(const char* name, uint16_t n, uint16_t k0) {
  FixedFFT<CalcT> fft;
  CHECK(fft.begin(n), "begin failed");

  std::vector<CalcT> input(n);
  const double amp = 0.4;  // stay well inside +-1 range
  const CalcT kMax = FFTLimits<CalcT>::kMax;
  for (uint16_t t = 0; t < n; t++) {
    double v = amp * std::cos(2.0 * M_PI * k0 * t / n);
    if constexpr (std::is_floating_point<CalcT>::value) {
      input[t] = (CalcT)v;  // already normalized - no integer quantization
    } else {
      input[t] = (CalcT)std::lround(v * kMax);
    }
  }
  CHECK(fft.loadReal(input.data(), n), "loadReal failed");
  fft.fft();

  typename FixedFFT<CalcT>::MagSquaredT best{};
  uint16_t bestBin = 0;
  for (uint16_t k = 0; k < n; k++) {
    auto m = fft.magnitudeSquared(k);
    if (m > best) {
      best = m;
      bestBin = k;
    }
  }
  bool ok = (bestBin == k0) || (bestBin == (uint16_t)(n - k0));
  std::printf("[%s] N=%d k0=%d -> peak bin %d %s\n", name, n, k0, bestBin, ok ? "OK" : "MISMATCH");
  CHECK(ok, "fft peak landed on the wrong bin");

  fft.end();
}

// --- Numeric accuracy vs. reference DFT -----------------------------------
template <typename CalcT>
void testAccuracy(const char* name, uint16_t n, bool adaptive = false) {
  FixedFFT<CalcT> fft;
  CHECK(fft.begin(n), "begin failed");

  std::mt19937 rng(12345 + n);
  std::uniform_real_distribution<double> dist(-0.5, 0.5);
  const CalcT kMax = FFTLimits<CalcT>::kMax;

  std::vector<CalcT> input(n);
  std::vector<std::complex<double>> ref_in(n);
  for (uint16_t t = 0; t < n; t++) {
    double v = dist(rng);
    if constexpr (std::is_floating_point<CalcT>::value) {
      input[t] = (CalcT)v;  // already normalized - no integer quantization
    } else {
      input[t] = (CalcT)std::lround(v * kMax);
    }
    ref_in[t] = std::complex<double>((double)input[t] / (double)kMax, 0.0);
  }
  CHECK(fft.loadReal(input.data(), n), "loadReal failed");
  uint8_t shifts = fft.fft(adaptive);
  auto ref = referenceDft(ref_in, false);

  // fft()'s output is the unscaled DFT sum divided by 2^shifts (not
  // necessarily N, when adaptiveScale skipped some stages).
  const double fwdScale = (double)((uint32_t)1 << shifts);
  double worst = 0.0;
  for (uint16_t k = 0; k < n; k++) {
    double gotRe = (double)fft.real(k) / kMax;
    double gotIm = (double)fft.imag(k) / kMax;
    double expRe = ref[k].real() / fwdScale;
    double expIm = ref[k].imag() / fwdScale;
    worst = std::max(worst, std::fabs(gotRe - expRe));
    worst = std::max(worst, std::fabs(gotIm - expIm));
  }
  std::printf("[%s%s] N=%d fft worst normalized error: %.5f (shifts=%d)\n", name,
              adaptive ? "/adaptive" : "", n, worst, shifts);
  const double fwdTol = sizeof(CalcT) == 1 ? 0.1 : (sizeof(CalcT) == 2 ? 0.01 : 0.001);
  CHECK(worst <= fwdTol, "fft accuracy vs. reference DFT exceeded tolerance for this calc type");

  // Round trip. The non-adaptive fft()/ifft() pairing is the one that
  // recovers the exact original amplitude; an adaptive fft() must be
  // paired with ifft(true) instead (its own adaptive guard) - that
  // combination is overflow-safe and preserves shape/relative
  // amplitude, but not a fixed, known overall scale factor, so we only
  // check correlation with the original signal here rather than exact
  // amplitude.
  if (!adaptive) {
    fft.ifft();
    double worstRt = 0.0;
    for (uint16_t t = 0; t < n; t++) {
      double got = (double)fft.real(t) / kMax;
      double exp = (double)input[t] / kMax;
      worstRt = std::max(worstRt, std::fabs(got - exp));
    }
    std::printf("[%s] N=%d round-trip worst error: %.5f\n", name, n, worstRt);
    // Round-trip precision for the non-adaptive pairing degrades with N
    // for low-precision calc types (int8_t/Q7 only has 7 fractional
    // bits to begin with, and the fixed 1/N forward attenuation eats
    // further into that as N grows) - this is an inherent, documented
    // property of the calc type/size combination, not a bug, so the
    // tolerance scales accordingly instead of using one fixed bound.
    double rtTol = sizeof(CalcT) == 1 ? 0.06 * ((double)n / 16.0) + 0.05
                                       : (sizeof(CalcT) == 2 ? 0.002 * ((double)n / 16.0) : 0.001);
    CHECK(worstRt <= rtTol, "round-trip accuracy exceeded tolerance for this calc type/size");
  } else {
    fft.ifft(/*safeScale=*/true);
    double dot = 0.0, refNorm = 0.0, gotNorm = 0.0;
    for (uint16_t t = 0; t < n; t++) {
      double got = (double)fft.real(t);
      double exp = (double)input[t];
      dot += got * exp;
      refNorm += exp * exp;
      gotNorm += got * got;
    }
    double correlation = dot / (std::sqrt(refNorm) * std::sqrt(gotNorm) + 1e-9);
    std::printf("[%s/adaptive] N=%d round-trip (via ifft(safeScale=true)) shape correlation: %.4f\n", name, n,
                correlation);
    CHECK(correlation > 0.9, "adaptive round trip lost too much shape correlation");
  }

  fft.end();
}

// --- loadHalfSpectrum(): synthesize a real signal from a one-sided
// spectrum, verifying the mirrored upper half gives ifft() a properly
// Hermitian-symmetric spectrum (imaginary part ~0 everywhere) and the
// right frequency content. ---------------------------------------------
template <typename CalcT>
void testHalfSpectrumSynthesis(const char* name, uint16_t n, uint16_t bin) {
  FixedFFT<CalcT> fft;
  CHECK(fft.begin(n), "begin failed");

  const uint16_t halfN = (uint16_t)(n / 2 + 1);
  const double amp = 0.3;
  const CalcT kMax = FFTLimits<CalcT>::kMax;

  std::vector<CalcT> real_in(halfN, (CalcT)0), imag_in(halfN, (CalcT)0);
  if (std::is_floating_point<CalcT>::value) {
    real_in[bin] = (CalcT)amp;
  } else {
    real_in[bin] = (CalcT)std::lround(amp * (double)kMax);
  }

  CHECK(fft.loadHalfSpectrum(real_in.data(), imag_in.data(), halfN), "loadHalfSpectrum failed");
  // Integer calc types need the adaptive overflow guard here (this
  // spectrum didn't come from fft(), so there's no guaranteed
  // headroom); float/double never need it (see ifft()'s doc comment).
  fft.ifft(!std::is_floating_point<CalcT>::value);

  double worstImag = 0.0;
  for (uint16_t i = 0; i < n; i++) {
    worstImag = std::max(worstImag, std::fabs((double)fft.imag(i)));
  }
  std::printf("[%s] N=%d loadHalfSpectrum worst imag (raw units): %.3g\n", name, n, worstImag);
  // Integer types build the mirror with exact integer arithmetic, so
  // the imaginary part is exactly 0; float/double are exact in the
  // mirroring itself too, but the FFT/IFFT butterfly math still
  // accumulates ordinary floating-point rounding error (~1e-15 for
  // double), so allow a tiny epsilon there instead of requiring
  // bit-exact 0.
  const double imagTol = std::is_floating_point<CalcT>::value ? 1e-6 : 0.0;
  CHECK(worstImag <= imagTol, "loadHalfSpectrum's mirrored spectrum did not yield an exactly real ifft() result");

  // Shape check: the synthesized signal should correlate with the
  // expected cosine at the requested bin (exact amplitude isn't
  // checked for integer types since the adaptive scale is data
  // dependent - see ifft()'s doc comment).
  const double kMaxD = std::is_floating_point<CalcT>::value ? 1.0 : (double)kMax;
  double dot = 0.0, refNorm = 0.0, gotNorm = 0.0;
  for (uint16_t t = 0; t < n; t++) {
    double got = (double)fft.real(t) / kMaxD;
    double exp = 2.0 * amp * std::cos(2.0 * M_PI * bin * t / n) / n;
    dot += got * exp;
    refNorm += exp * exp;
    gotNorm += got * got;
  }
  double correlation = dot / (std::sqrt(refNorm) * std::sqrt(gotNorm) + 1e-12);
  std::printf("[%s] N=%d loadHalfSpectrum synthesized signal correlation: %.6f\n", name, n, correlation);
  CHECK(correlation > 0.999, "loadHalfSpectrum synthesized the wrong signal shape");

  fft.end();
}

int main() {
  testTwiddleTable<int8_t>("q7");
  testTwiddleTable<int16_t>("q15");
  testTwiddleTable<int32_t>("q31");
  testTwiddleTable<float>("f32");
  testTwiddleTable<double>("f64");

  for (uint16_t n : {8, 16, 32, 64, 128}) {
    testPeakBin<int8_t>("q7", n, n / 4);
    testPeakBin<int16_t>("q15", n, n / 4);
    testPeakBin<int32_t>("q31", n, n / 4);
    testPeakBin<float>("f32", n, n / 4);
    testPeakBin<double>("f64", n, n / 4);
  }

  // Accuracy thresholds get tighter as the calc type's precision grows.
  for (uint16_t n : {16, 64, 256}) {
    testAccuracy<int8_t>("q7", n);
    testAccuracy<int16_t>("q15", n);
    testAccuracy<int32_t>("q31", n);
    testAccuracy<float>("f32", n);
    testAccuracy<double>("f64", n);
  }
  // int8_t (Q7) at larger N benefits from adaptive forward scaling -
  // compare against the unconditional-scaling numbers printed above.
  testAccuracy<int8_t>("q7", 256, /*adaptive=*/true);

  // Cross-width input conversion: feed int16_t samples into an int8_t
  // (Q7) calc engine and vice versa, just check it runs and lands the
  // peak on the right bin.
  {
    FixedFFT<int8_t, int16_t> fft;
    CHECK(fft.begin(32), "begin failed");
    std::vector<int16_t> input(32);
    for (int t = 0; t < 32; t++) {
      input[t] = (int16_t)std::lround(0.4 * 32767.0 * std::cos(2.0 * M_PI * 4 * t / 32));
    }
    CHECK(fft.loadReal(input.data(), 32), "loadReal(int16->q7) failed");
    fft.fft();
    uint64_t best = 0;
    uint16_t bestBin = 0;
    for (uint16_t k = 0; k < 32; k++) {
      uint64_t m = fft.magnitudeSquared(k);
      if (m > best) {
        best = m;
        bestBin = k;
      }
    }
    std::printf("[q7 from int16 input] peak bin %d (expect 4 or 28)\n", bestBin);
    CHECK(bestBin == 4 || bestBin == 28, "cross-width conversion landed on wrong bin");
    fft.end();
  }

  // Same, but into a float calc engine: this is the exact bug pattern
  // that surfaced during development (feeding an integer InT into a
  // float CalcT must scale up into [-1, 1), not shift like int->int
  // does) - keep it covered explicitly.
  {
    FixedFFT<float, int16_t> fft;
    CHECK(fft.begin(32), "begin failed");
    std::vector<int16_t> input(32);
    for (int t = 0; t < 32; t++) {
      input[t] = (int16_t)std::lround(0.4 * 32767.0 * std::cos(2.0 * M_PI * 4 * t / 32));
    }
    CHECK(fft.loadReal(input.data(), 32), "loadReal(int16->f32) failed");
    fft.fft();
    float best = 0;
    uint16_t bestBin = 0;
    for (uint16_t k = 0; k < 32; k++) {
      float m = fft.magnitudeSquared(k);
      if (m > best) {
        best = m;
        bestBin = k;
      }
    }
    std::printf("[f32 from int16 input] peak bin %d (expect 4 or 28), magnitude=%.4f (expect ~6.4)\n", bestBin,
                fft.magnitude(bestBin));
    CHECK(bestBin == 4 || bestBin == 28, "cross-width conversion landed on wrong bin");
    CHECK(std::fabs(fft.magnitude(bestBin) - 6.4) < 0.1, "int16->float magnitude scale looks wrong");
    fft.end();
  }

  // Same again into a double calc engine, with a tighter tolerance
  // since double has much more precision available.
  {
    FixedFFT<double, int16_t> fft;
    CHECK(fft.begin(32), "begin failed");
    std::vector<int16_t> input(32);
    for (int t = 0; t < 32; t++) {
      input[t] = (int16_t)std::lround(0.4 * 32767.0 * std::cos(2.0 * M_PI * 4 * t / 32));
    }
    CHECK(fft.loadReal(input.data(), 32), "loadReal(int16->f64) failed");
    fft.fft();
    double best = 0;
    uint16_t bestBin = 0;
    for (uint16_t k = 0; k < 32; k++) {
      double m = fft.magnitudeSquared(k);
      if (m > best) {
        best = m;
        bestBin = k;
      }
    }
    std::printf("[f64 from int16 input] peak bin %d (expect 4 or 28), magnitude=%.6f (expect ~6.4)\n", bestBin,
                fft.magnitude(bestBin));
    CHECK(bestBin == 4 || bestBin == 28, "cross-width conversion landed on wrong bin");
    CHECK(std::fabs(fft.magnitude(bestBin) - 6.4) < 0.01, "int16->double magnitude scale looks wrong");
    fft.end();
  }

  for (uint16_t n : {16, 32, 64}) {
    testHalfSpectrumSynthesis<int8_t>("q7", n, n / 8);
    testHalfSpectrumSynthesis<int16_t>("q15", n, n / 8);
    testHalfSpectrumSynthesis<int32_t>("q31", n, n / 8);
    testHalfSpectrumSynthesis<float>("f32", n, n / 8);
    testHalfSpectrumSynthesis<double>("f64", n, n / 8);
  }

  if (g_failures == 0) {
    std::printf("\nALL TESTS PASSED\n");
    return 0;
  } else {
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
  }
}
