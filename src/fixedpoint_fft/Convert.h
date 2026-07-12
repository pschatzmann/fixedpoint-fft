#pragma once
#include <stdint.h>

#include "Traits.h"

namespace fixedpoint_fft {

/// Rescales a sample from an arbitrary input type into the representation
/// used by CalcT. InT is normally int8_t/int16_t/int32_t (a full-range
/// signed PCM sample), but a floating-point InT is also accepted - it is
/// assumed to already be normalized to this library's [-1, 1) full-scale
/// convention (the same convention loadReal()/load() produce for a float
/// CalcT), so it needs scaling only if CalcT itself is an integer type.
///
/// Between two integer types, this shifts into/out of the extra headroom
/// bits to land in CalcT's Qn range. Between an integer InT and a float
/// CalcT (or vice versa), it normalizes by/multiplies by 2^bits instead
/// of shifting. This only ever runs at load time (O(N), once per
/// transform), never inside the butterfly hot loop, so using plain
/// int64_t/double intermediates here costs nothing that matters while
/// staying trivially overflow-safe for every combination of input/calc
/// type.
///
/// Implemented via a class template partial-specialized on the
/// InT/CalcT "is floating point" flags (rather than "if constexpr",
/// which needs C++17) so it stays buildable under the AVR toolchain's
/// plain C++11 - see the comment at the top of FixedFFTEngine.h.
template <typename InT, typename CalcT, bool kInFloat, bool kOutFloat>
struct FFTConvertHelper;

/// FFTConvertHelper specialization: int -> int, shift into/out of the
/// extra headroom bits.
template <typename InT, typename CalcT>
struct FFTConvertHelper<InT, CalcT, false, false> {
  static inline CalcT convert(InT v) {
    const int kInBits = (int)sizeof(InT) * 8 - 1;
    const int kOutBits = (int)sizeof(CalcT) * 8 - 1;
    if (kOutBits >= kInBits) {
      return (CalcT)(((int64_t)v) << (kOutBits - kInBits));
    } else {
      return (CalcT)(((int64_t)v) >> (kInBits - kOutBits));
    }
  }
};

/// FFTConvertHelper specialization: int -> float, normalize into [-1, 1).
template <typename InT, typename CalcT>
struct FFTConvertHelper<InT, CalcT, false, true> {
  static inline CalcT convert(InT v) {
    const int kInBits = (int)sizeof(InT) * 8 - 1;
    const CalcT kInvScale = (CalcT)1 / (CalcT)((int64_t)1 << kInBits);
    return (CalcT)v * kInvScale;
  }
};

/// FFTConvertHelper specialization: float -> int, scale a
/// [-1, 1)-normalized value up into CalcT's Qn range.
template <typename InT, typename CalcT>
struct FFTConvertHelper<InT, CalcT, true, false> {
  static inline CalcT convert(InT v) {
    const int kOutBits = (int)sizeof(CalcT) * 8 - 1;
    const double kScale = (double)(((uint64_t)1 << kOutBits) - 1);
    return (CalcT)((double)v * kScale);
  }
};

/// FFTConvertHelper specialization: float -> float, already in the same
/// convention, pass through.
template <typename InT, typename CalcT>
struct FFTConvertHelper<InT, CalcT, true, true> {
  static inline CalcT convert(InT v) { return (CalcT)v; }
};

/// Converts a single sample of type InT into CalcT's representation -
/// see the FFTConvertHelper doc comment above for the conversion rules.
template <typename InT, typename CalcT>
inline CalcT fftConvertSample(InT v) {
  return FFTConvertHelper<InT, CalcT, FFTIsFloat<InT>::value, FFTIsFloat<CalcT>::value>::convert(v);
}

}  // namespace fixedpoint_fft
