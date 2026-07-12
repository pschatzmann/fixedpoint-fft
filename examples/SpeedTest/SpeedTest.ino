// Benchmarks fft()/ifft() speed for every supported calculation type -
// int8_t, int16_t, int32_t, float, and (off AVR) double - so you can
// compare the actual speed/precision tradeoff on your specific board
// instead of taking it on faith.
//
// kN defaults to 64 so this fits comfortably in an AVR Uno's 2KB of
// RAM even for the largest calc type (int32_t/double: 64 * 4 * 2 =
// 512 bytes). On a board with more RAM/headroom (ESP32, SAMD, ...),
// bump kN up for a more meaningful comparison - at N=64 the fastest
// calc types on fast MCUs may finish in just a few microseconds,
// close to micros()'s own resolution.

#include "FixedFFT.h"

constexpr uint16_t kN = 64;            // FFT size used for the benchmark
constexpr uint16_t kIterations = 200;  // repetitions averaged per measurement

// Shared input signal (content doesn't matter for timing) - int16_t so
// it can feed every calc type via loadReal()'s automatic rescaling.
int16_t g_samples[kN];

template <typename CalcT>
void benchmark(const char* name) {
  FixedFFT<CalcT, int16_t> fft;  // g_samples below is int16_t regardless of CalcT
  if (!fft.begin(kN)) {
    Serial.print(name);
    Serial.println("\tbegin() failed (kN exceeds this calc type's kMaxN?)");
    return;
  }

  unsigned long fftTotalMicros = 0;
  unsigned long ifftTotalMicros = 0;

  for (uint16_t iter = 0; iter < kIterations; iter++) {
    fft.loadReal(g_samples, kN);

    unsigned long t0 = micros();
    fft.fft();
    fftTotalMicros += micros() - t0;

    unsigned long t1 = micros();
    fft.ifft();
    ifftTotalMicros += micros() - t1;
  }

  float fftAvg = (float)fftTotalMicros / kIterations;
  float ifftAvg = (float)ifftTotalMicros / kIterations;

  Serial.print(name);
  Serial.print("\t");
  Serial.print((int)sizeof(CalcT));
  Serial.print("\t");
  Serial.print(fftAvg, 2);
  Serial.print("\t");
  Serial.print(ifftAvg, 2);
  Serial.print("\t");
  Serial.println(1.0e6f / fftAvg, 0);

  fft.end();
}

void setup() {
  Serial.begin(115200);
  while (!Serial) {}

  for (uint16_t t = 0; t < kN; t++) {
    g_samples[t] = (int16_t)(10000.0f * cosf(2.0f * PI * 5 * t / kN));
  }

  Serial.print("Benchmarking N=");
  Serial.print(kN);
  Serial.print(", ");
  Serial.print(kIterations);
  Serial.println(" iterations per calc type");
  Serial.println("type\tbytes\tfft()us\tifft()us\tfft()/s");

  benchmark<int8_t>("int8_t");
  benchmark<int16_t>("int16_t");
  benchmark<int32_t>("int32_t");
  benchmark<float>("float");
  // double is a desktop-only precision option (see README) - the AVR
  // toolchain doesn't even support FixedFFT<double> (its `double` is
  // bit-identical to `float`, so there's nothing genuine to benchmark),
  // but it's worth including on every other target.
#if !defined(ARDUINO_ARCH_AVR)
  benchmark<double>("double");
#endif

  Serial.println("done");
}

void loop() {}
