#pragma once
#include <stdint.h>
#include <stddef.h>
#include <math.h>

#include "Convert.h"
#include "Traits.h"
#include "Twiddles.h"

// This header must stay buildable with plain C++11: the AVR toolchain
// Arduino ships (avr-gcc, via -std=gnu++11) has neither a C++17
// compiler nor a C++ standard library beyond the bare language (no
// <type_traits>, no <cmath>, not even a "c++" include directory) - only
// the C library's <math.h>, which is used below for the float/double
// calc types' magnitude(). So instead of "if constexpr" (C++17) to
// branch FixedFFTEngineCore's behavior between the integer Qn calc
// types and float/double, this file uses plain C++11 explicit
// specialization: the class body below is the generic (integer)
// implementation, and fft()/ifft()/transform()/magnitude()/
// magnitudeSquared() are each re-specialized for CalcT=float and
// CalcT=double immediately after the class, with their own, separately
// written bodies.
//
// FixedFFTEngineCore<CalcT> holds everything that depends only on the
// calculation type. FixedFFTEngine<CalcT, InT> (the public API, aliased
// as FixedFFT) derives from it and adds the input-sample loaders
// (loadReal/load/loadHalfSpectrum), which depend on a second,
// independent type: the type your raw samples are stored as. Splitting
// it this way means the float/double specializations below only have
// to be written once per CalcT, rather than once per (CalcT, InT)
// combination.

namespace fixedpoint_fft {

/// In-place, iterative radix-2 decimation-in-time fixed-point FFT core -
/// everything that depends only on the calculation type. See
/// FixedFFTEngine (below) for the public, two-template-parameter API
/// that adds the input-sample loaders on top of this.
///
///   CalcT - the type all math is performed in: int8_t (Q7, default -
///           fastest / least RAM on 8-bit AVR), int16_t (Q15), int32_t
///           (Q31), float, or double.
///
///           float is worth choosing specifically on MCUs with a
///           hardware single-precision FPU (ESP32/ESP32-S3, ARM
///           Cortex-M4F/M7, ...): a plain float multiply-add is already
///           correctly scaled (no Qn rescale-by-shift step, unlike the
///           integer types), floats have enough dynamic range that the
///           per-stage overflow-guard scaling is unnecessary, and on
///           that hardware it's about as fast as int32_t while being
///           easier to reason about. On integer-only targets like plain
///           AVR, float is emulated in software and is much slower than
///           any of the integer calc types - stick to int8_t/int16_t/
///           int32_t there.
///
///           double is a desktop-only precision option: real MCU targets
///           either have no hardware double (AVR, ESP32 - it's emulated
///           in software, no faster than float there, sometimes slower)
///           or, on ARM Cortex-M4F/M7, no benefit over float (the FPU is
///           single-precision only). On an actual desktop/server x86-64
///           CPU, though, double is native and roughly as fast as float
///           while giving ~15-17 significant decimal digits instead of
///           ~7 - useful for e.g. host-side reference/verification work
///           (this library's own test suite uses it that way) or
///           precision-critical offline analysis.
///
/// Storage: exactly two CalcT arrays of length N (real + imaginary) - no
/// other buffer is allocated by the transform itself. You can supply your
/// own static/stack arrays via begin(n, real, imag) for a genuinely
/// zero-heap, in-place transform, or let begin(n) allocate them for you.
///
/// No sin()/cos() from the C/C++ library is ever called: twiddle factors
/// come from the offline-generated flash tables in generated/. magnitude()
/// uses an integer bit-by-bit square root for the integer calc types; for
/// CalcT=float/double it uses sqrtf()/sqrt() (<math.h>), which on an
/// FPU-equipped target is a single hardware instruction - faster than
/// reimplementing the integer algorithm for a type it wasn't designed for.
template <typename CalcT = int8_t>
class FixedFFTEngineCore {
 public:
  using Wide = typename FFTWide<CalcT>::type;
  using MagSquaredT = typename FFTConditional<FFTIsFloat<CalcT>::value, CalcT, uint64_t>::type;
  static constexpr uint8_t kQBits = FFTQBits<CalcT>::value;
  static constexpr uint16_t kMaxN = FFTTwiddleTraits<CalcT>::kMaxN;

  FixedFFTEngineCore() = default;
  ~FixedFFTEngineCore() { end(); }

  /// Use externally supplied, already-allocated buffers (true in-place,
  /// zero extra heap/stack use by the library). Both arrays must have at
  /// least n elements and n must be a power of two <= kMaxN.
  bool begin(uint16_t n, CalcT* real_buf, CalcT* imag_buf) {
    if (!isPow2(n) || n < 2 || n > kMaxN || real_buf == nullptr || imag_buf == nullptr) {
      return false;
    }
    end();
    len = n;
    log2len = ilog2(n);
    p_real = real_buf;
    p_imag = imag_buf;
    owns_buffers = false;
    return true;
  }

  /// Convenience overload: allocates the two working buffers on the heap.
  /// Costs 2 * n * sizeof(CalcT) bytes, freed again by end() / the
  /// destructor.
  bool begin(uint16_t n) {
    if (!isPow2(n) || n < 2 || n > kMaxN) return false;
    // Arduino cores (AVR, ESP32, SAMD, ...) implement global operator
    // new on top of malloc() and return nullptr on failure rather than
    // throwing, so a plain new[] here is already effectively "nothrow".
    CalcT* r = new CalcT[n];
    CalcT* i = new CalcT[n];
    if (r == nullptr || i == nullptr) {
      delete[] r;
      delete[] i;
      return false;
    }
    if (!begin(n, r, i)) {
      delete[] r;
      delete[] i;
      return false;
    }
    owns_buffers = true;
    return true;
  }

  /// Releases any heap-allocated buffers (no-op for externally supplied
  /// buffers) and resets this instance so begin() can be called again.
  void end() {
    if (owns_buffers) {
      delete[] p_real;
      delete[] p_imag;
    }
    p_real = nullptr;
    p_imag = nullptr;
    len = 0;
    log2len = 0;
    owns_buffers = false;
  }

  /// Returns the FFT length this instance was configured with via begin().
  uint16_t size() const { return len; }

  /// Forward FFT, fully in place. By default every one of the log2(N)
  /// stages unconditionally divides its butterfly outputs by 2, so
  /// intermediate values can never overflow regardless of input
  /// amplitude; the result is X[k] * (1/N) in CalcT's Qn format. This is
  /// the fastest option and is exact enough for int16_t/int32_t calc
  /// types at any supported N.
  ///
  /// For int8_t (Q7) at larger N, unconditionally halving every stage
  /// can waste precision on stages whose values never needed it, leaving
  /// too few bits for a usable round trip. Pass adaptiveScale=true to
  /// only halve a stage when a buffer scan shows it's actually needed
  /// (extra O(N) scan per stage, better precision, same-or-less overall
  /// attenuation).
  ///
  /// Returns the total number of halvings actually applied (always
  /// log2(N) when adaptiveScale=false): fft()'s output equals the
  /// standard unscaled DFT sum divided by 2^(return value), so
  /// left-shifting the result by that amount recovers the unscaled
  /// spectrum.
  ///
  /// ifft()'s default (fast) mode assumes the paired fft() call used
  /// exactly log2(N) halvings, i.e. adaptiveScale=false - that pairing
  /// is the one to use when you need the exact original amplitude back.
  /// If you pass adaptiveScale=true here, pair it with ifft(true)
  /// (its own adaptive guard) rather than the default ifft(): the
  /// combination is overflow-safe and preserves the signal's shape /
  /// relative amplitude, but the two adaptive passes don't cancel each
  /// other's scaling exactly, so the overall output amplitude is
  /// data-dependent rather than a fixed, known factor. Do not try to
  /// compensate by shifting the ifft() output yourself - by the time
  /// ifft() returns, any precision an insufficient shift count lost
  /// to saturation mid-transform is already gone.
  ///
  /// CalcT=float/double each have their own specialization of this
  /// method (see below the class) with different behavior: adaptiveScale
  /// is ignored and it always returns 0, since per-stage scaling is
  /// never needed for a floating-point calc type.
  uint8_t fft(bool adaptiveScale = false) {
    uint8_t shifts = 0;
    transform(false, /*forceScale=*/!adaptiveScale, adaptiveScale, &shifts);
    if (!adaptiveScale) shifts = log2len;
    return shifts;
  }

  /// Inverse FFT, fully in place. By default this does *not* rescale per
  /// stage: that exactly undoes fft()'s 1/N attenuation for the common
  /// fft()->ifft() round trip (recovers the original signal, up to
  /// fixed-point rounding) and is the fastest option.
  ///
  /// If you need to run ifft() on an arbitrary/full-scale spectrum that
  /// didn't come from this library's fft() - so headroom isn't
  /// guaranteed - pass safeScale=true to enable an adaptive per-stage
  /// overflow guard (O(N) buffer scan per stage, only shifts when a
  /// stage's values actually need it). Returns the number of halvings
  /// actually applied (0 when safeScale=false, or when no stage needed
  /// one); multiply the result by 2^(return value) to compensate.
  ///
  /// CalcT=float/double each have their own specialization of this
  /// method (see below the class): safeScale is ignored, and it always
  /// performs one exact, single-pass division by N at the end (cheaper
  /// and more accurate than per-stage halving), since fft() never
  /// attenuates for a floating-point calc type.
  uint8_t ifft(bool safeScale = false) {
    uint8_t shifts = 0;
    transform(true, /*forceScale=*/false, /*adaptiveScale=*/safeScale, &shifts);
    return shifts;
  }

  /// Returns the real part of bin i in CalcT's units.
  CalcT real(uint16_t i) const { return p_real[i]; }
  /// Returns the imaginary part of bin i in CalcT's units.
  CalcT imag(uint16_t i) const { return p_imag[i]; }
  /// Returns a pointer to the real-part working buffer.
  CalcT* realData() { return p_real; }
  /// Returns a pointer to the imaginary-part working buffer.
  CalcT* imagData() { return p_imag; }

  /// Squared magnitude in raw Qn units (i.e. re^2 + im^2, computed
  /// without narrowing) - avoids a sqrt() call entirely when only
  /// relative bin comparisons are needed (e.g. finding the peak bin).
  /// CalcT=float/double each have their own specialization returning a
  /// plain float/double.
  MagSquaredT magnitudeSquared(uint16_t i) const {
    int64_t re = p_real[i];
    int64_t im = p_imag[i];
    return (uint64_t)(re * re + im * im);
  }

  /// Magnitude in the same Qn units as real()/imag(), via an integer
  /// bit-by-bit square root (no <cmath>, no sqrt()). CalcT=float has its
  /// own specialization using sqrtf()/sqrt() instead (see below the
  /// class).
  CalcT magnitude(uint16_t i) const {
    return fftSaturate<CalcT>((Wide)fftIsqrt(magnitudeSquared(i)));
  }

 private:
  CalcT* p_real = nullptr;
  CalcT* p_imag = nullptr;
  uint16_t len = 0;
  uint8_t log2len = 0;
  bool owns_buffers = false;

  static bool isPow2(uint16_t n) { return n != 0 && (n & (n - 1)) == 0; }
  static uint8_t ilog2(uint16_t n) {
    uint8_t r = 0;
    while (n > 1) {
      n >>= 1;
      r++;
    }
    return r;
  }

  /// Standard iterative in-place bit-reversal permutation.
  void bitReverse() {
    uint16_t j = 0;
    for (uint16_t i = 1; i < len; i++) {
      uint16_t bit = len >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) {
        CalcT tr = p_real[i];
        p_real[i] = p_real[j];
        p_real[j] = tr;
        CalcT ti = p_imag[i];
        p_imag[i] = p_imag[j];
        p_imag[j] = ti;
      }
    }
  }

  /// True if any current sample would overflow CalcT once summed with
  /// another same-magnitude sample, i.e. it's already using more than
  /// half the type's range. Used only by ifft(safeScale=true).
  bool needsScale() const {
    const Wide thresh = (Wide)FFTLimits<CalcT>::kMax >> 1;
    for (uint16_t i = 0; i < len; i++) {
      Wide ar = p_real[i];
      if (ar < 0) ar = -ar;
      Wide ai = p_imag[i];
      if (ai < 0) ai = -ai;
      if (ar > thresh || ai > thresh) return true;
    }
    return false;
  }

  /// CalcT=float/double each have their own specialization of this
  /// method (see below the class) with a simpler, unscaled butterfly
  /// core.
  void transform(bool inverse, bool forceScale, bool adaptiveScale, uint8_t* shiftsOut) {
    bitReverse();
    const uint16_t N = len;

    for (uint16_t size = 2; size <= N; size <<= 1) {
      const uint16_t half = size >> 1;
      const uint16_t angleStep = kMaxN / size;

      bool scaleThisStage = forceScale;
      if (!forceScale && adaptiveScale) {
        scaleThisStage = needsScale();
        if (scaleThisStage && shiftsOut != nullptr) (*shiftsOut)++;
      }

      for (uint16_t base = 0; base < N; base += size) {
        uint16_t angleIdx = 0;
        for (uint16_t k = 0; k < half; k++, angleIdx = (uint16_t)(angleIdx + angleStep)) {
          CalcT c, s;
          fftTwiddle<CalcT>(angleIdx, c, s);
          const CalcT tw_re = c;
          const CalcT tw_im = inverse ? s : (CalcT)(-s);

          const uint16_t i1 = base + k;
          const uint16_t i2 = i1 + half;

          const CalcT a_re = p_real[i2];
          const CalcT a_im = p_imag[i2];

          const Wide ac = (Wide)a_re * tw_re;
          const Wide ad = (Wide)a_re * tw_im;
          const Wide bc = (Wide)a_im * tw_re;
          const Wide bd = (Wide)a_im * tw_im;

          // Q2n complex product, rescaled back down to Qn (round to
          // nearest).
          const Wide half_r = (Wide)1 << (kQBits - 1);
          const Wide t_re = (ac - bd + half_r) >> kQBits;
          const Wide t_im = (ad + bc + half_r) >> kQBits;

          const Wide sum_re = (Wide)p_real[i1] + t_re;
          const Wide sum_im = (Wide)p_imag[i1] + t_im;
          const Wide diff_re = (Wide)p_real[i1] - t_re;
          const Wide diff_im = (Wide)p_imag[i1] - t_im;

          if (scaleThisStage) {
            p_real[i1] = fftRoundShift<CalcT>(sum_re, 1);
            p_imag[i1] = fftRoundShift<CalcT>(sum_im, 1);
            p_real[i2] = fftRoundShift<CalcT>(diff_re, 1);
            p_imag[i2] = fftRoundShift<CalcT>(diff_im, 1);
          } else {
            p_real[i1] = fftSaturate<CalcT>(sum_re);
            p_imag[i1] = fftSaturate<CalcT>(sum_im);
            p_real[i2] = fftSaturate<CalcT>(diff_re);
            p_imag[i2] = fftSaturate<CalcT>(diff_im);
          }
        }
      }
    }
  }
};

// ---------------------------------------------------------------------
// CalcT=float specializations. Written as separate, self-contained
// bodies (C++11 explicit member specialization) rather than branching
// inside the generic methods above, since the AVR toolchain has no
// "if constexpr" (C++17) available - see the file-level comment.
// ---------------------------------------------------------------------

/// FixedFFTEngineCore<float>::transform specialization: a plain,
/// unscaled butterfly core (no Qn rescale-by-shift, no overflow guard).
template <>
inline void FixedFFTEngineCore<float>::transform(bool inverse, bool /*forceScale*/, bool /*adaptiveScale*/,
                                                  uint8_t* /*shiftsOut*/) {
  bitReverse();
  const uint16_t N = len;

  for (uint16_t size = 2; size <= N; size <<= 1) {
    const uint16_t half = size >> 1;
    const uint16_t angleStep = kMaxN / size;

    for (uint16_t base = 0; base < N; base += size) {
      uint16_t angleIdx = 0;
      for (uint16_t k = 0; k < half; k++, angleIdx = (uint16_t)(angleIdx + angleStep)) {
        float c, s;
        fftTwiddle<float>(angleIdx, c, s);
        const float tw_re = c;
        const float tw_im = inverse ? s : -s;

        const uint16_t i1 = base + k;
        const uint16_t i2 = i1 + half;

        const float a_re = p_real[i2];
        const float a_im = p_imag[i2];

        // A plain float multiply is already correctly scaled - no Qn
        // rescale-by-shift step, no overflow guard, needed.
        const float t_re = a_re * tw_re - a_im * tw_im;
        const float t_im = a_re * tw_im + a_im * tw_re;

        const float sum_re = p_real[i1] + t_re;
        const float sum_im = p_imag[i1] + t_im;
        const float diff_re = p_real[i1] - t_re;
        const float diff_im = p_imag[i1] - t_im;

        p_real[i1] = sum_re;
        p_imag[i1] = sum_im;
        p_real[i2] = diff_re;
        p_imag[i2] = diff_im;
      }
    }
  }
}

/// FixedFFTEngineCore<float>::fft specialization: adaptiveScale is
/// ignored - per-stage scaling is never needed for float. Output is the
/// plain unscaled DFT sum (return value 0), the convention most other
/// FFT libraries (FFTW, numpy, ESP-DSP, ...) use for float.
template <>
inline uint8_t FixedFFTEngineCore<float>::fft(bool /*adaptiveScale*/) {
  transform(false, false, false, nullptr);
  return 0;
}

/// FixedFFTEngineCore<float>::ifft specialization: safeScale is ignored -
/// since fft() never attenuates for float, this always performs one
/// exact, single-pass division by N at the end (cheaper and more
/// accurate than per-stage halving) so that fft()->ifft() recovers the
/// original signal directly.
template <>
inline uint8_t FixedFFTEngineCore<float>::ifft(bool /*safeScale*/) {
  transform(true, false, false, nullptr);
  const float invN = 1.0f / (float)len;
  for (uint16_t i = 0; i < len; i++) {
    p_real[i] = p_real[i] * invN;
    p_imag[i] = p_imag[i] * invN;
  }
  return log2len;
}

/// FixedFFTEngineCore<float>::magnitudeSquared specialization: plain
/// float re^2 + im^2, no integer narrowing.
template <>
inline float FixedFFTEngineCore<float>::magnitudeSquared(uint16_t i) const {
  return p_real[i] * p_real[i] + p_imag[i] * p_imag[i];
}

/// FixedFFTEngineCore<float>::magnitude specialization: sqrtf()
/// (<math.h>) is a single hardware instruction on an FPU-equipped
/// target - faster than repurposing the integer bit-by-bit square root
/// for a type it wasn't designed for.
template <>
inline float FixedFFTEngineCore<float>::magnitude(uint16_t i) const {
  return sqrtf(magnitudeSquared(i));
}

// ---------------------------------------------------------------------
// CalcT=double specializations - a straight mirror of the float ones
// above, for the desktop-only precision option described in the class
// doc comment. Not selected by anything in this library automatically;
// only used if you explicitly write FixedFFT<double>.
//
// Excluded on AVR: unlike an ordinary template, an explicit
// specialization's body is compiled as soon as the compiler sees it,
// whether or not anything ever calls it - so transform()'s call to
// fftTwiddle<double>() would require fftTableRead(const double*, ...)
// to exist even in a sketch that never uses CalcT=double. AVR
// deliberately has no such overload (see Traits.h - avr-gcc's `double`
// is bit-identical to `float` there anyway, so there's nothing genuine
// to support), which would otherwise break every AVR build regardless
// of which calc type it actually uses.
// ---------------------------------------------------------------------
#if !defined(ARDUINO_ARCH_AVR)

/// FixedFFTEngineCore<double>::transform specialization: a plain,
/// unscaled butterfly core (no Qn rescale-by-shift, no overflow guard).
template <>
inline void FixedFFTEngineCore<double>::transform(bool inverse, bool /*forceScale*/, bool /*adaptiveScale*/,
                                                   uint8_t* /*shiftsOut*/) {
  bitReverse();
  const uint16_t N = len;

  for (uint16_t size = 2; size <= N; size <<= 1) {
    const uint16_t half = size >> 1;
    const uint16_t angleStep = kMaxN / size;

    for (uint16_t base = 0; base < N; base += size) {
      uint16_t angleIdx = 0;
      for (uint16_t k = 0; k < half; k++, angleIdx = (uint16_t)(angleIdx + angleStep)) {
        double c, s;
        fftTwiddle<double>(angleIdx, c, s);
        const double tw_re = c;
        const double tw_im = inverse ? s : -s;

        const uint16_t i1 = base + k;
        const uint16_t i2 = i1 + half;

        const double a_re = p_real[i2];
        const double a_im = p_imag[i2];

        // A plain double multiply is already correctly scaled - no Qn
        // rescale-by-shift step, no overflow guard, needed.
        const double t_re = a_re * tw_re - a_im * tw_im;
        const double t_im = a_re * tw_im + a_im * tw_re;

        const double sum_re = p_real[i1] + t_re;
        const double sum_im = p_imag[i1] + t_im;
        const double diff_re = p_real[i1] - t_re;
        const double diff_im = p_imag[i1] - t_im;

        p_real[i1] = sum_re;
        p_imag[i1] = sum_im;
        p_real[i2] = diff_re;
        p_imag[i2] = diff_im;
      }
    }
  }
}

/// FixedFFTEngineCore<double>::fft specialization: adaptiveScale is
/// ignored - per-stage scaling is never needed for double. Output is the
/// plain unscaled DFT sum (return value 0), the convention most other
/// FFT libraries (FFTW, numpy, ESP-DSP, ...) use for floating-point calc
/// types.
template <>
inline uint8_t FixedFFTEngineCore<double>::fft(bool /*adaptiveScale*/) {
  transform(false, false, false, nullptr);
  return 0;
}

/// FixedFFTEngineCore<double>::ifft specialization: safeScale is ignored -
/// since fft() never attenuates for double, this always performs one
/// exact, single-pass division by N at the end (cheaper and more
/// accurate than per-stage halving) so that fft()->ifft() recovers the
/// original signal directly.
template <>
inline uint8_t FixedFFTEngineCore<double>::ifft(bool /*safeScale*/) {
  transform(true, false, false, nullptr);
  const double invN = 1.0 / (double)len;
  for (uint16_t i = 0; i < len; i++) {
    p_real[i] = p_real[i] * invN;
    p_imag[i] = p_imag[i] * invN;
  }
  return log2len;
}

/// FixedFFTEngineCore<double>::magnitudeSquared specialization: plain
/// double re^2 + im^2, no integer narrowing.
template <>
inline double FixedFFTEngineCore<double>::magnitudeSquared(uint16_t i) const {
  return p_real[i] * p_real[i] + p_imag[i] * p_imag[i];
}

/// FixedFFTEngineCore<double>::magnitude specialization: sqrt()
/// (<math.h>) is a single hardware instruction on a target with native
/// double (e.g. desktop x86-64) - faster than repurposing the integer
/// bit-by-bit square root for a type it wasn't designed for.
template <>
inline double FixedFFTEngineCore<double>::magnitude(uint16_t i) const {
  return sqrt(magnitudeSquared(i));
}

#endif  // !defined(ARDUINO_ARCH_AVR)

/// The public FFT engine: FixedFFTEngineCore<CalcT> (everything that
/// depends only on the calculation type) plus the input-sample loaders,
/// which depend on a second, independent type - InT, the type your raw
/// samples are stored as. InT defaults to CalcT, so existing code that
/// only ever names one template argument (e.g. FixedFFT<int8_t>) keeps
/// working unchanged, loading samples already in that same type.
///
/// Naming a different InT lets you load samples in one representation
/// while computing in another without any per-call casting or template
/// argument on the load call itself - e.g. FixedFFT<int8_t, int16_t>
/// computes at int8_t (Q7) speed/RAM while accepting raw int16_t ADC
/// samples directly:
///
///   FixedFFT<int8_t, int16_t> fft;
///   fft.begin(64);
///   int16_t samples[64] = { /* ... */ };
///   fft.loadReal(samples, 64);   // rescaled from int16_t into Q7 on load
///
/// InT may be int8_t, int16_t, int32_t, float, or double, independent of
/// CalcT - see fftConvertSample (Convert.h) for the conversion rules
/// between every combination.
template <typename CalcT = int8_t, typename InT = CalcT>
class FixedFFTEngine : public FixedFFTEngineCore<CalcT> {
 public:
  /// Loads a real-valued signal (imaginary part cleared to 0), rescaling
  /// from InT's full range into CalcT's Qn range. InT is normally
  /// int8_t/int16_t/int32_t; a floating-point InT is also accepted and
  /// is assumed to already be normalized to this library's [-1, 1)
  /// convention (see fftConvertSample).
  bool loadReal(const InT* input, uint16_t n) {
    if (n != this->size() || input == nullptr) return false;
    CalcT* re = this->realData();
    CalcT* im = this->imagData();
    for (uint16_t i = 0; i < n; i++) {
      re[i] = fftConvertSample<InT, CalcT>(input[i]);
      im[i] = 0;
    }
    return true;
  }

  /// Loads a complex-valued signal.
  bool load(const InT* real_in, const InT* imag_in, uint16_t n) {
    if (n != this->size() || real_in == nullptr || imag_in == nullptr) return false;
    CalcT* re = this->realData();
    CalcT* im = this->imagData();
    for (uint16_t i = 0; i < n; i++) {
      re[i] = fftConvertSample<InT, CalcT>(real_in[i]);
      im[i] = fftConvertSample<InT, CalcT>(imag_in[i]);
    }
    return true;
  }

  /// Loads a one-sided spectrum - bins 0..N/2 (halfN = N/2 + 1 complex
  /// values) - and mirrors it into the redundant upper half (bins
  /// N/2+1..N-1) via conjugate symmetry (real[N-k] = real[k],
  /// imag[N-k] = -imag[k]) before a subsequent ifft() runs.
  ///
  /// This is the standard way to synthesize a real time-domain signal
  /// from a frequency-domain specification: ifft() only produces a
  /// purely real result (imaginary part ~0 at every bin, up to rounding)
  /// when the spectrum has this Hermitian symmetry. Building the full
  /// N-length spectrum by hand and calling load() works too, but this
  /// halves the data you need to specify and guarantees the symmetry is
  /// exact rather than something you have to get right yourself.
  ///
  /// Bins 0 and N/2 (Nyquist - always present since N is a power of two,
  /// so always even) map to themselves under conjugate symmetry, so
  /// their imaginary part is forced to exactly 0 regardless of what you
  /// pass in - a nonzero value there would otherwise make ifft()'s
  /// result subtly complex.
  bool loadHalfSpectrum(const InT* real_in, const InT* imag_in, uint16_t halfN) {
    const uint16_t n = this->size();
    if (halfN != (uint16_t)(n / 2 + 1) || real_in == nullptr || imag_in == nullptr) return false;
    CalcT* re = this->realData();
    CalcT* im = this->imagData();
    for (uint16_t k = 0; k < halfN; k++) {
      re[k] = fftConvertSample<InT, CalcT>(real_in[k]);
      im[k] = fftConvertSample<InT, CalcT>(imag_in[k]);
    }
    im[0] = 0;
    im[n / 2] = 0;
    for (uint16_t k = 1; k < n / 2; k++) {
      re[n - k] = re[k];
      im[n - k] = fftSaturate<CalcT>((typename FFTWide<CalcT>::type) - (typename FFTWide<CalcT>::type)im[k]);
    }
    return true;
  }
};

}  // namespace fixedpoint_fft
