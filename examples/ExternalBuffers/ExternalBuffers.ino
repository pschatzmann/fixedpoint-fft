// Zero-heap, fully in-place FFT: the real/imaginary working buffers are
// plain static arrays supplied to begin(), so the library itself never
// calls new/malloc. Uses int16_t (Q15) as the calculation type for
// better precision than the int8_t default, fed from int16_t samples
// (e.g. a typical 16-bit ADC/I2S reading buffer) named via FixedFFT's
// second template parameter.

#include "FixedFFT.h"

constexpr uint16_t kN = 128;
constexpr uint32_t kSampleRate = 16000;  // Hz - the rate sensor_samples was captured at

FixedFFT<int16_t, int16_t> fft;  // int16_t (Q15) calc type, int16_t input samples
static int16_t real_buf[kN];
static int16_t imag_buf[kN];
static int16_t sensor_samples[kN];  // e.g. a 16-bit ADC reading buffer

// Bin k covers the frequency k * kSampleRate / kN.
float binFrequency(uint16_t k) { return (float)k * kSampleRate / kN; }

void setup() {
  Serial.begin(115200);
  while (!Serial) {}

  if (!fft.begin(kN, real_buf, imag_buf)) {
    Serial.println("fft.begin() failed");
    return;
  }

  for (uint16_t t = 0; t < kN; t++) {
    float v = 0.3f * cosf(2.0f * PI * 12 * t / kN) + 0.2f * cosf(2.0f * PI * 30 * t / kN);
    sensor_samples[t] = (int16_t)lroundf(v * 32767.0f);
  }

  fft.loadReal(sensor_samples, kN);
  fft.fft();

  uint16_t peakBin = 0;
  uint64_t peakMag = 0;
  for (uint16_t k = 0; k < kN / 2; k++) {
    uint64_t m = fft.magnitudeSquared(k);
    if (m > peakMag) {
      peakMag = m;
      peakBin = k;
    }
  }
  Serial.print("peak bin: ");
  Serial.print(peakBin);
  Serial.print("\tfrequency(Hz): ");
  Serial.println(binFrequency(peakBin));
}

void loop() {}
