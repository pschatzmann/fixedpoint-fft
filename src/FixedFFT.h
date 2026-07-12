#pragma once

#include "fixedpoint_fft/FixedFFTEngine.h"

namespace fixedpoint_fft {

/// fixedpoint-fft: a fast, in-place FFT for Arduino, built primarily
/// around fixed-point integer math, with float/double also available as
/// calc types on hardware where they make sense.
///
///   - No sin()/cos() calls, ever, for any calc type. Twiddle factors
///     come from a small offline-generated flash-resident table (see
///     tools/generate_twiddles.py). magnitude() uses an integer square
///     root for the integer calc types, or sqrtf()/sqrt() (<math.h>)
///     for float/double.
///   - Fully in place: exactly two CalcT arrays of length N (real +
///     imaginary), which you may supply yourself for zero heap use.
///   - Two independent template parameters: FixedFFT<CalcT, InT>. CalcT
///     (default int8_t, ~7 bits of precision, fastest/least RAM) is the
///     type all math is performed in - int8_t, int16_t, int32_t, float,
///     or double. InT (defaults to CalcT) is the type your raw samples
///     are stored as - also any of those five, independent of CalcT -
///     so e.g. FixedFFT<int8_t, int16_t> computes at int8_t speed while
///     accepting raw int16_t samples directly, no per-call casting.
///
/// Basic usage:
///
///   #include "FixedFFT.h"
///   fixedpoint_fft::FixedFFT<int8_t> fft;   // int8_t calc + input type (default)
///
///   fft.begin(64);                          // allocate 2*64 bytes
///   fft.loadReal(samples, 64);               // samples: int8_t[64]
///   fft.fft();
///   for (int i = 0; i < 32; i++) {
///     int8_t mag = fft.magnitude(i);
///   }
///
/// See examples/ for complete sketches, and README.md for the math
/// behind the fixed-point scaling and how to regenerate the twiddle
/// tables for larger FFT sizes.
template <typename CalcT = int8_t, typename InT = CalcT>
using FixedFFT = FixedFFTEngine<CalcT, InT>;

}  // namespace fixedpoint_fft

// Arduino sketches are compiled as if everything were in one global
// translation unit, so pull the namespace in automatically there -
// sketches can just write FixedFFT<int8_t> instead of
// fixedpoint_fft::FixedFFT<int8_t>. Regular C++ projects that #include
// this header outside of Arduino keep the namespace explicit.
#ifdef ARDUINO
using namespace fixedpoint_fft;
#endif
