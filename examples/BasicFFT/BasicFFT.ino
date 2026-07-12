// Basic forward FFT + magnitude spectrum, using the default int8_t
// (Q7) calculation type - fastest, smallest RAM footprint - fed from
// int16_t samples (e.g. a typical 16-bit ADC/I2S reading buffer),
// named via FixedFFT's second template parameter.
//
// Generates a synthetic tone in setup(), runs the FFT once, and prints
// the magnitude of each bin over Serial.

#include "FixedFFT.h"

constexpr uint16_t kN = 64;
constexpr uint32_t kSampleRate = 8000;  // Hz - the rate `samples` was captured at
constexpr uint16_t kToneBin = 8;        // bin the synthetic tone should land on

FixedFFT<int8_t, int16_t> fft;  // int8_t calc type, int16_t input samples
int16_t samples[kN];

// Bin k covers the frequency k * kSampleRate / kN.
float binFrequency(uint16_t k) { return (float)k * kSampleRate / kN; }

void setup() {
  Serial.begin(115200);
  while (!Serial) {}

  if (!fft.begin(kN)) {
    Serial.println("fft.begin() failed - kN must be a power of two <= FixedFFT<int8_t>::kMaxN");
    return;
  }

  // A single-frequency test tone at ~40% of full scale.
  for (uint16_t t = 0; t < kN; t++) {
    float v = 0.4f * cosf(2.0f * PI * kToneBin * t / kN);
    samples[t] = (int16_t)lroundf(v * 32767.0f);
  }

  fft.loadReal(samples, kN);
  fft.fft();

  Serial.println("bin\tfrequency(Hz)\tmagnitude");
  for (uint16_t k = 0; k < kN / 2; k++) {
    Serial.print(k);
    Serial.print('\t');
    Serial.print(binFrequency(k));
    Serial.print('\t');
    Serial.println(fft.magnitude(k));
  }
}

void loop() {}
