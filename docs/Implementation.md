# Implementation details

Background on how fixedpoint-fft works internally - not needed to just
*use* the library (see [README.md](README.md) for that), but useful if
you're curious how the "no trig calls, no floating point required"
promise is actually kept, or if you're modifying the library yourself.

## How it avoids trigonometric calls

`tools/generate_twiddles.py` is a host-side Python script (never
compiled into a sketch) that computes `cos(2*pi*i / max_n)` for one
quarter turn (`i = 0 .. max_n/4`), quantizes it into the target Qn
fixed-point format (or leaves it as a plain float/double for the
`float`/`double` calc types - no quantization needed there), and writes
it out as a `const`/`PROGMEM` C array under `src/fixedpoint_fft/generated/`
(one table per calc type: `twiddles_q7.h`, `twiddles_q15.h`,
`twiddles_q31.h`, `twiddles_f32.h`, `twiddles_f64.h`). At runtime,
`Twiddles.h` reconstructs the full circle
from that single quarter table via quadrant folding (mirroring/negating,
no interpolation, no recomputation), and one table (sized for `kMaxN`,
4096 by default) serves every power-of-two FFT length up to that
maximum - smaller transforms just use a larger stride into the same
table.

If you need FFT sizes larger than the checked-in default supports,
rerun the generator with a bigger `--max-n` (must be a power of two)
and it will overwrite the generated headers in place:

```sh
python3 tools/generate_twiddles.py --max-n 4096
```

## Fixed-point scaling model

Each `CalcT` is treated as a signed Qn fraction in `[-1, 1)`:
`int8_t` -> Q7, `int16_t` -> Q15, `int32_t` -> Q31. All multiply-
accumulate steps happen in a doubled-width accumulator (`int16_t` for
`int8_t`, `int32_t` for `int16_t`, `int64_t` for `int32_t` - the same
sizing CMSIS-DSP uses for its q7/q15/q31 kernels), so a single butterfly
can never overflow; only the final store back into `CalcT` narrows
(with saturation as a last-resort safety net, never expected to trigger
in normal operation).

- **`fft()`** divides every stage's butterfly outputs by 2
  unconditionally. This guarantees no overflow regardless of input
  amplitude and is the fastest option; the result is the standard DFT
  sum scaled by `1/N`. It returns the number of halvings actually
  applied (`log2(N)`) so you can left-shift back to the unscaled
  spectrum if you need it.
- **`ifft()`** by default does **not** rescale per stage - that exactly
  undoes `fft()`'s `1/N` attenuation for the common `fft()` -> `ifft()`
  round trip (recovers the original signal, up to fixed-point rounding)
  and is the fastest option. This assumes the spectrum came from this
  library's `fft()` (or otherwise has headroom); running `ifft()`
  directly on an arbitrary full-scale spectrum can overflow.
- Both `fft(true)` and `ifft(true)` switch to **adaptive** scaling: an
  O(N) buffer scan before each stage decides whether that stage
  actually needs halving, applying it only when necessary. This trades
  some speed for precision - useful for `int8_t` at larger N, where
  unconditional halving can throw away more bits than it needs to. Only
  pair an adaptive `fft(true)` with an adaptive `ifft(true)` - mixing
  adaptive forward scaling with the default (non-adaptive) inverse will
  reconstruct the signal at an unpredictable, likely-clipped amplitude
  (see the doc comment on `fft()` in `FixedFFTEngine.h` for why).

## Why float/double needed their own code path

The AVR toolchain Arduino ships (`avr-gcc`) compiles C++ as
`-std=gnu++11` and has no C++ standard library beyond the bare
language - no `<type_traits>`, no `<cmath>`, not even a `c++` include
directory - only the C library's `<math.h>`. That rules out the
obvious C++17 way to branch this library's behavior between the
integer Qn calc types and `float` (`if constexpr` + `std::is_floating_point`),
since `int8_t`/`int16_t`/`int32_t` (the types AVR sketches actually
use) still need to compile under that same C++11 restriction.

Instead, `FixedFFTEngineCore`'s `transform()`/`fft()`/`ifft()`/
`magnitude()`/`magnitudeSquared()` are written once as the generic
(integer) implementation, and re-specialized for `CalcT=float`/
`CalcT=double` via plain C++11 explicit member-function specialization
immediately below the class in `FixedFFTEngine.h` - each with its own,
separately written body, rather than one shared body with a
runtime/constexpr branch inside it. Two small hand-rolled helpers
(`FFTIsFloat`, `FFTConditional` in `Traits.h`) stand in for
`std::is_floating_point`/`std::conditional` where a compile-time type
decision is still needed (e.g. `fftConvertSample`'s InT/CalcT
combinations in `Convert.h`, or the `magnitudeSquared()` return type),
without depending on `<type_traits>`.

One more wrinkle specific to `double`: an *explicit* specialization's
body is compiled as soon as the compiler sees it, whether or not
anything ever calls it - unlike an ordinary template, which is only
instantiated on use. So the `CalcT=double` specializations are wrapped
in `#if !defined(ARDUINO_ARCH_AVR)` in `FixedFFTEngine.h`: without that
guard, `transform()`'s call to `fftTwiddle<double>()` would require a
`double` PROGMEM table reader to exist on AVR even in a sketch that
never uses `CalcT=double`, breaking every AVR build regardless of which
calc type it actually uses.

This is also why `CalcT` and `InT` (the calc type and the input sample
type) are split across two classes rather than living as two
parameters on one: `FixedFFTEngineCore<CalcT>` holds everything above
(the numeric engine, keyed only on `CalcT`, including the float/double
specializations); the public `FixedFFTEngine<CalcT, InT>` (aliased as
`FixedFFT`) derives from it and adds `loadReal()`/`load()`/
`loadHalfSpectrum()`, which are the only things that actually depend on
`InT`. Explicit specialization is keyed on the *exact* template
argument list, so if `transform()` etc. lived directly on a
two-parameter class, the float/double specializations would need one
copy per `(CalcT, InT)` combination instead of one per `CalcT` -
`FixedFFTEngineCore<CalcT>` avoids that combinatorial blow-up entirely.
