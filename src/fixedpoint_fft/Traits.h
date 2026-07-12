#pragma once
#include <stdint.h>

/// PROGMEM / flash-table read abstraction. On AVR, const arrays default to
/// RAM unless explicitly marked PROGMEM and read back with pgm_read_*; on
/// every other Arduino core (ESP32, SAMD, RP2040, Teensy, ...) plain const
/// arrays already live in flash and are addressable directly, so the
/// attribute and the accessor both become no-ops there.
#if defined(ARDUINO_ARCH_AVR)
#include <avr/pgmspace.h>
#define FFT_TABLE_ATTR PROGMEM
#else
#define FFT_TABLE_ATTR
#endif

namespace fixedpoint_fft {

/// Hand-rolled equivalents of std::is_floating_point/std::conditional.
/// The AVR toolchain Arduino ships (avr-gcc) has no C++ standard library
/// at all - no <type_traits>, no <cmath>, nothing beyond the bare
/// language - so this library cannot depend on <type_traits> even though
/// every other target (ESP32, ARM, ...) has it. These are the only two
/// pieces of it actually needed, and both are trivial to write directly.
template <typename T>
struct FFTIsFloat {
  static constexpr bool value = false;
};
/// FFTIsFloat specialization: float is a floating point type.
template <>
struct FFTIsFloat<float> {
  static constexpr bool value = true;
};
/// FFTIsFloat specialization: double is a floating point type.
template <>
struct FFTIsFloat<double> {
  static constexpr bool value = true;
};

/// Hand-rolled equivalent of std::conditional - see FFTIsFloat above for
/// why this can't just be <type_traits>.
template <bool B, typename T, typename F>
struct FFTConditional {
  using type = T;
};
/// FFTConditional specialization for the false branch.
template <typename T, typename F>
struct FFTConditional<false, T, F> {
  using type = F;
};

#if defined(ARDUINO_ARCH_AVR)
/// Reads an int8_t twiddle table entry out of AVR PROGMEM (flash).
inline int8_t fftTableRead(const int8_t* table, uint16_t idx) {
  return (int8_t)pgm_read_byte(&table[idx]);
}
/// Reads an int16_t twiddle table entry out of AVR PROGMEM (flash).
inline int16_t fftTableRead(const int16_t* table, uint16_t idx) {
  return (int16_t)pgm_read_word(&table[idx]);
}
/// Reads an int32_t twiddle table entry out of AVR PROGMEM (flash).
inline int32_t fftTableRead(const int32_t* table, uint16_t idx) {
  return (int32_t)pgm_read_dword(&table[idx]);
}
/// Reads a float twiddle table entry out of AVR PROGMEM (flash).
inline float fftTableRead(const float* table, uint16_t idx) { return pgm_read_float(&table[idx]); }
/// No PROGMEM double reader: avr-gcc's `double` is bit-identical to
/// `float` (no real double-precision hardware or format on AVR), so
/// CalcT=double isn't meaningfully supported there anyway - use
/// CalcT=float on AVR instead. Trying to instantiate FixedFFT<double> on
/// AVR fails to compile here with a clear "no matching function" error
/// rather than silently reinterpreting bytes.
#else
/// Reads an int8_t twiddle table entry (plain flash-resident array).
inline int8_t fftTableRead(const int8_t* table, uint16_t idx) { return table[idx]; }
/// Reads an int16_t twiddle table entry (plain flash-resident array).
inline int16_t fftTableRead(const int16_t* table, uint16_t idx) { return table[idx]; }
/// Reads an int32_t twiddle table entry (plain flash-resident array).
inline int32_t fftTableRead(const int32_t* table, uint16_t idx) { return table[idx]; }
/// Reads a float twiddle table entry (plain flash-resident array).
inline float fftTableRead(const float* table, uint16_t idx) { return table[idx]; }
/// Reads a double twiddle table entry (plain flash-resident array).
inline double fftTableRead(const double* table, uint16_t idx) { return table[idx]; }
#endif

/// Doubles the bit width of a calc type so multiply-accumulate steps never
/// overflow: int8 x int8 -> int16, int16 x int16 -> int32, int32 x int32 ->
/// int64. This mirrors how CMSIS-DSP sizes its q7/q15/q31 accumulators.
/// float maps to itself: a plain single-precision multiply is already
/// correctly scaled (no Qn renormalization step is needed the way fixed
/// point needs one), and float has enough dynamic range that widening
/// buys nothing.
template <typename T>
struct FFTWide;
/// FFTWide specialization: int8_t widens to int16_t.
template <>
struct FFTWide<int8_t> {
  using type = int16_t;
};
/// FFTWide specialization: int16_t widens to int32_t.
template <>
struct FFTWide<int16_t> {
  using type = int32_t;
};
/// FFTWide specialization: int32_t widens to int64_t.
template <>
struct FFTWide<int32_t> {
  using type = int64_t;
};
/// FFTWide specialization: float widens to itself (no widening needed).
template <>
struct FFTWide<float> {
  using type = float;
};
/// FFTWide specialization: double widens to itself (no widening needed).
template <>
struct FFTWide<double> {
  using type = double;
};

/// Number of fractional bits in the Qn representation of a calc type, i.e.
/// int8_t -> Q7, int16_t -> Q15, int32_t -> Q31 (sign bit is not counted).
/// Meaningless for float (there's no Qn rescale step for it - see
/// FixedFFTEngine::transform()), defined only so the value compiles.
template <typename T>
struct FFTQBits;
/// FFTQBits specialization: int8_t is Q7.
template <>
struct FFTQBits<int8_t> {
  static constexpr uint8_t value = 7;
};
/// FFTQBits specialization: int16_t is Q15.
template <>
struct FFTQBits<int16_t> {
  static constexpr uint8_t value = 15;
};
/// FFTQBits specialization: int32_t is Q31.
template <>
struct FFTQBits<int32_t> {
  static constexpr uint8_t value = 31;
};
/// FFTQBits specialization: unused/meaningless for float.
template <>
struct FFTQBits<float> {
  static constexpr uint8_t value = 0;
};
/// FFTQBits specialization: unused/meaningless for double.
template <>
struct FFTQBits<double> {
  static constexpr uint8_t value = 0;
};

/// Largest/smallest representable Qn value for an integer calc type.
template <typename T>
struct FFTLimits {
  static constexpr T kMax = (T)((((typename FFTWide<T>::type)1) << FFTQBits<T>::value) - 1);
  static constexpr T kMin = -kMax - 1;
};

/// float/double use the same nominal [-1, 1) full-scale convention as
/// the integer Qn types by convention (fftConvertSample() normalizes
/// into it), but nothing here enforces or clamps to it - see
/// fftSaturate<float>/fftSaturate<double> below. These constants exist
/// purely so generic code that reads FFTLimits<CalcT>::kMax/kMin
/// (including this library's own tests) works uniformly across every
/// CalcT.
template <>
struct FFTLimits<float> {
  static constexpr float kMax = 1.0f;
  static constexpr float kMin = -1.0f;
};
/// FFTLimits specialization for double - see the float specialization
/// above for why these values exist without being enforced.
template <>
struct FFTLimits<double> {
  static constexpr double kMax = 1.0;
  static constexpr double kMin = -1.0;
};

/// Saturating narrow of a wide accumulator back down to CalcT.
template <typename CalcT>
inline CalcT fftSaturate(typename FFTWide<CalcT>::type v) {
  if (v > FFTLimits<CalcT>::kMax) return FFTLimits<CalcT>::kMax;
  if (v < FFTLimits<CalcT>::kMin) return FFTLimits<CalcT>::kMin;
  return (CalcT)v;
}

/// float/double have enormous dynamic range compared to the handful of
/// bits an int8_t/int16_t/int32_t Qn value gets, so there's no realistic
/// overflow to guard against here - narrowing would only risk
/// incorrectly clipping a legitimately large value (e.g. an unscaled
/// forward transform's DC bin, which can be much larger than 1.0). No
/// clamping.
template <>
inline float fftSaturate<float>(float v) {
  return v;
}
/// fftSaturate specialization for double - see the float specialization
/// above for the rationale.
template <>
inline double fftSaturate<double>(double v) {
  return v;
}

/// Arithmetic right shift with round-to-nearest, then saturate to CalcT.
/// Used to rescale a Q(2n) multiply result back to Qn, and to apply the
/// per-stage 1/2 butterfly scaling that keeps the in-place FFT overflow-free.
template <typename CalcT>
inline CalcT fftRoundShift(typename FFTWide<CalcT>::type v, uint8_t shift) {
  if (shift == 0) return fftSaturate<CalcT>(v);
  typename FFTWide<CalcT>::type half = (typename FFTWide<CalcT>::type)1 << (shift - 1);
  return fftSaturate<CalcT>((v + half) >> shift);
}

/// Integer square root (bit-by-bit, non-restoring) - used for magnitude()
/// so the library never has to pull in <cmath>/sqrt() either. Takes a
/// 64-bit input so it is safe to reuse for every CalcT's magnitude-squared
/// range (up to ~2^63 for the int32_t/Q31 case).
inline uint32_t fftIsqrt(uint64_t v) {
  uint64_t res = 0;
  uint64_t bit = (uint64_t)1 << 62;
  while (bit > v) bit >>= 2;
  while (bit != 0) {
    if (v >= res + bit) {
      v -= res + bit;
      res = (res >> 1) + bit;
    } else {
      res >>= 1;
    }
    bit >>= 2;
  }
  return (uint32_t)res;
}

}  // namespace fixedpoint_fft
