#pragma once
#include <stdint.h>

#include "Traits.h"
#include "generated/twiddles_f32.h"
#include "generated/twiddles_f64.h"
#include "generated/twiddles_q7.h"
#include "generated/twiddles_q15.h"
#include "generated/twiddles_q31.h"

/// Reconstructs full-circle twiddle factors from a single quarter-wave
/// (0 .. pi/2) cosine table per calc type, generated offline by
/// tools/generate_twiddles.py. No sin()/cos()/any trig library call ever
/// happens here or anywhere else in this library - it is pure integer
/// table lookup plus quadrant sign/mirror logic.
///
/// One table (sized for kTwiddle*MaxN) serves every power-of-two FFT
/// length up to that maximum: smaller transforms simply index the table
/// with a larger stride.

namespace fixedpoint_fft {

/// Per-calc-type twiddle table accessor: exposes the table's maximum
/// supported FFT length (kMaxN) and a pointer to its quarter-wave data.
template <typename CalcT>
struct FFTTwiddleTraits;

/// FFTTwiddleTraits specialization for int8_t (Q7).
template <>
struct FFTTwiddleTraits<int8_t> {
  static constexpr uint16_t kMaxN = kTwiddleQ7MaxN;
  static const int8_t* table() { return kTwiddleQ7Quarter; }
};

/// FFTTwiddleTraits specialization for int16_t (Q15).
template <>
struct FFTTwiddleTraits<int16_t> {
  static constexpr uint16_t kMaxN = kTwiddleQ15MaxN;
  static const int16_t* table() { return kTwiddleQ15Quarter; }
};

/// FFTTwiddleTraits specialization for int32_t (Q31).
template <>
struct FFTTwiddleTraits<int32_t> {
  static constexpr uint16_t kMaxN = kTwiddleQ31MaxN;
  static const int32_t* table() { return kTwiddleQ31Quarter; }
};

/// FFTTwiddleTraits specialization for float.
template <>
struct FFTTwiddleTraits<float> {
  static constexpr uint16_t kMaxN = kTwiddleF32MaxN;
  static const float* table() { return kTwiddleF32Quarter; }
};

/// FFTTwiddleTraits specialization for double.
template <>
struct FFTTwiddleTraits<double> {
  static constexpr uint16_t kMaxN = kTwiddleF64MaxN;
  static const double* table() { return kTwiddleF64Quarter; }
};

/// Returns cos(2*pi*idx/kMaxN) and sin(2*pi*idx/kMaxN) for idx in [0, kMaxN).
/// The caller (FixedFFTEngine) combines these into the actual rotation
/// factor and picks the sign of the imaginary part for forward vs. inverse.
template <typename CalcT>
inline void fftTwiddle(uint16_t idx, CalcT& cosOut, CalcT& sinOut) {
  using Traits = FFTTwiddleTraits<CalcT>;
  constexpr uint16_t Q = Traits::kMaxN / 4;
  const CalcT* t = Traits::table();

  uint16_t quadrant = idx / Q;
  uint16_t r = idx % Q;

  switch (quadrant) {
    case 0:
      cosOut = fftTableRead(t, r);
      sinOut = fftTableRead(t, Q - r);
      break;
    case 1:
      cosOut = (CalcT)(-fftTableRead(t, Q - r));
      sinOut = fftTableRead(t, r);
      break;
    case 2:
      cosOut = (CalcT)(-fftTableRead(t, r));
      sinOut = (CalcT)(-fftTableRead(t, Q - r));
      break;
    default:  // quadrant 3
      cosOut = fftTableRead(t, Q - r);
      sinOut = (CalcT)(-fftTableRead(t, r));
      break;
  }
}

}  // namespace fixedpoint_fft
