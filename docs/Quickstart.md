## Quick start

### Forward FFT

`FixedFFT<CalcT, InT>` takes two template parameters: `CalcT` is the
type all math is performed in, `InT` (defaults to `CalcT`) is the type
your raw samples are stored as. Naming them independently - e.g.
`int8_t` calc type fed from `int16_t` ADC/I2S samples, as below - avoids
any per-call casting or conversion of your own: `loadReal()` rescales
from `InT` into `CalcT`'s Qn range on the way in.

```cpp
#include "FixedFFT.h"

fixedpoint_fft::FixedFFT<int8_t, int16_t> fft;  // int8_t calc type, int16_t input samples

void setup() {
  fft.begin(64);                 // allocates 2 * 64 bytes internally
}

void loop() {
  int16_t samples[64] = { /* ... */ };
  fft.loadReal(samples, 64);     // rescaled from int16_t into Q7; imag cleared to 0
  fft.fft();
  for (int k = 0; k < 32; k++) {
    int8_t mag = fft.magnitude(k);  // no sqrt() call - integer isqrt
  }
}
```

(If your samples are already in `CalcT`'s own representation, just
omit the second template argument - `FixedFFT<int8_t>` is equivalent to
`FixedFFT<int8_t, int8_t>`.)

For a genuinely zero-heap, fully in-place transform, supply your own
buffers instead of letting `begin(n)` allocate them (these are always
`CalcT`-typed working buffers, regardless of `InT`):

```cpp
static int8_t real_buf[64], imag_buf[64];
fft.begin(64, real_buf, imag_buf);
```

### Inverse FFT

You don't need to call `fft()` first to use `ifft()` - e.g. to
synthesize a time-domain signal from a spectrum you built yourself (a
handful of bins set to a magnitude/phase, everything else left at 0).
Two methods load a complex-valued spectrum (or signal) directly:

- **`load(real_in, imag_in, n)`** loads all `n` complex bins as given,
  no assumptions made. Use this when you want full control - e.g. the
  spectrum is genuinely complex-valued, or you're loading a
  time-domain signal that already has a nonzero imaginary part.
- **`loadHalfSpectrum(real_in, imag_in, halfN)`** loads only the
  non-redundant half of a spectrum - bins `0..N/2`, `halfN = N/2 + 1`
  complex values - and mirrors it into the upper half (bins
  `N/2+1..N-1`) via conjugate symmetry (`real[N-k] = real[k]`,
  `imag[N-k] = -imag[k]`) before `ifft()` runs. `ifft()` only produces a
  purely real result if the spectrum has this symmetry, so this is the
  method to use for the common case of synthesizing a real signal from
  a frequency-domain specification - you only specify `N/2 + 1` values
  instead of building the full symmetric spectrum by hand.

Since a spectrum loaded either way didn't come from this library's
`fft()` (so it has no guaranteed headroom), pass `safeScale=true` so
`ifft()` applies its adaptive per-stage overflow guard instead of
assuming one:

```cpp
int16_t real_in[33] = {0}, imag_in[33] = {0};  // N/2 + 1 = 33 for N=64
real_in[4] = 12000;                 // a single frequency component at bin 4
fft.loadHalfSpectrum(real_in, imag_in, 33);  // mirrors bins 33..63 automatically
fft.ifft(/*safeScale=*/true);       // inverse transform, in place
for (int i = 0; i < 64; i++) {
  int8_t sample = fft.real(i);      // synthesized time-domain signal
  // fft.imag(i) is ~0 at every sample, by construction
}
```

See `examples/SynthesizeFromSpectrum` for a complete sketch.

See `examples/` for complete sketches (`BasicFFT`, `ExternalBuffers`,
`InverseFFT`, `FloatFFT`, `SpeedTest`, `SynthesizeFromSpectrum`).

For how the library avoids trigonometric calls entirely and how its
fixed-point scaling model works, see
[IMPLEMENTATION.md](IMPLEMENTATION.md).
