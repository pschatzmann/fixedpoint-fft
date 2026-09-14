

## Choosing `float`/`double`

`float` is worth it specifically on FPU-equipped MCUs (ESP32/ESP32-S3,
ARM Cortex-M4F/M7, ...): a plain float multiply-add is already
correctly scaled, so it skips the rescale-by-shift step the integer
types need after every complex multiply, and it never needs the
saturating overflow guard applied to every butterfly output either -
float's huge dynamic range makes that overflow a non-issue, while the
integer types' narrow range makes it a real one. That's not a minor
saving: on hardware-FPU targets it consistently makes `float` several
times faster than `int32_t`, not merely comparable to it - e.g. on the
ESP32 (Xtensa LX6) an N=64 FFT takes 75.14 µs in `float` vs. 398.47 µs
in `int32_t`, and even the smaller `int16_t` calc type (250.49 µs) is
well behind. See [Performance.md](Performance.md) for the full
cross-device comparison. As a result, `fft()` always returns `0` (plain
unscaled DFT sum) and `magnitude()` uses `sqrtf()` instead of the
integer bit-by-bit square root for `float`. On plain AVR (no FPU),
`float` is emulated in software and much slower than the integer calc
types - stick to `int8_t`/`int16_t`/`int32_t` there. See
`examples/FloatFFT`.

`double` should be avoided on microcontrollers - it's unsupported on
AVR and gives no speed benefit over `float` elsewhere. It's meant for
desktop/server builds, where the extra precision (~15-17 significant
digits instead of `float`'s ~7) is useful for e.g. host-side reference/
verification work - this library's own test suite uses it that way.
