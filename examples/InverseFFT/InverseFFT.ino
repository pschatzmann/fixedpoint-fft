// Forward FFT followed by inverse FFT, demonstrating the round trip
// that recovers the original signal. Uses int32_t (Q31) as the
// calculation type since this sketch checks the reconstruction against
// the original samples and wants headroom to show a tight match, fed
// from int16_t samples (e.g. a typical 16-bit ADC/I2S reading buffer)
// named via FixedFFT's second template parameter.

#include "FixedFFT.h"

constexpr uint16_t kN = 32;
constexpr uint32_t kSampleRate = 4000;  // Hz - the rate samples was captured at

FixedFFT<int32_t, int16_t> fft;  // int32_t calc type, int16_t input samples
int16_t samples[kN];

// Bin k covers the frequency k * kSampleRate / kN.
float binFrequency(uint16_t k) { return (float)k * kSampleRate / kN; }

void setup() {
  Serial.begin(115200);
  while (!Serial) {}

  fft.begin(kN);

  for (uint16_t t = 0; t < kN; t++) {
    float v = 0.5f * cosf(2.0f * PI * 5 * t / kN);
    samples[t] = (int16_t)lroundf(v * 32767.0f);
  }

  fft.loadReal(samples, kN);
  fft.fft();  // forward: unconditionally scaled by 1/N, no overflow risk

  // Report the dominant frequency found before overwriting the buffer
  // with the inverse transform.
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

  fft.ifft();  // inverse: default (unscaled) mode exactly undoes fft()'s 1/N

  // samples[] is int16_t (Q15) and fft.real() is int32_t (Q31) - very
  // different raw ranges for the same normalized amplitude, so both are
  // printed normalized to [-1, 1) for a meaningful side-by-side compare.
  Serial.println("sample\toriginal\treconstructed");
  for (uint16_t t = 0; t < kN; t++) {
    Serial.print(t);
    Serial.print('\t');
    Serial.print(samples[t] / 32768.0f, 5);
    Serial.print('\t');
    Serial.println(fft.real(t) / 2147483648.0f, 5);
  }
}

void loop() {}
