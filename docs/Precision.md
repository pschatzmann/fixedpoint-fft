## Precision vs. calculation type - know the tradeoff

`int8_t` (Q7) gives you the smallest memory footprint and the fastest
per-butterfly math, but only ~7 bits of precision to start with. Its
accuracy holds up fine for peak/frequency detection (bin location) at
any supported N, but the **round-trip** `fft()`/`ifft()` precision
degrades as N grows, because the fixed `1/N` forward attenuation eats
further into an already-small budget of bits. As a rough guide from
this library's own test suite (`test/host_test.cpp`):

| CalcT     | N=16 | N=32 | N=64 | N=128 | N=256 | N=512 | N=1024 | N=2048 |
|-----------|------|------|------|-------|-------|-------|--------|--------|
| `int8_t`  | ~0.06 | ~0.17 | ~0.24 | ~0.53 | ~0.99 (unusable) | unusable | unusable | unusable |
| `int16_t` | ~0.0002 | ~0.0005 | ~0.001 | ~0.002 | ~0.004 | ~0.007 | ~0.016 | ~0.032 |
| `int32_t` | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 |
| `float`   | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 |
| `double`  | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 | ~0.000000 |

(Round-trip error, normalized to the `[-1, 1)` Qn range.) `int8_t`
becomes unusable for round trips from N=256 up - the fixed `1/N`
attenuation has eaten all 7 fractional bits by then, so the numbers
there are dominated by noise rather than a meaningful trend. If you
need `int8_t` speed/RAM at larger N, use `fft(true)`/`ifft(true)`
(adaptive scaling) rather than the defaults, or switch to `int16_t`.
Forward-only spectrum analysis (magnitude/bin detection, no `ifft()`)
is far less sensitive to this and works well with `int8_t` at any N.

For a detailed performance comparision by processor see [this table](Performance.md)
