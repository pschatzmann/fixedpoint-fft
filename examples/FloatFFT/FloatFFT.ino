// Forward FFT using the float (single-precision) calculation type.
//
// float is worth choosing specifically on MCUs with a hardware
// single-precision FPU - ESP32/ESP32-S3, ARM Cortex-M4F/M7 (e.g.
// Teensy 3.x/4.x, many STM32 "F4"/"F7" boards), RP2350's Cortex-M33 -
// where it runs about as fast as int32_t while avoiding the Qn
// rescale-by-shift step the integer calc types need, and gives far more
// headroom (no per-stage overflow scaling to think about). On an
// integer-only MCU like plain AVR, float is emulated in software and
// much slower than int8_t/int16_t/int32_t - use one of those there
// instead (see BasicFFT/ExternalBuffers).
//
// Twiddle factors still come from the offline-generated flash table
// (see tools/generate_twiddles.py) - no sin()/cos() call happens here
// either.

#include "FixedFFT.h"

constexpr uint16_t kN = 256;
constexpr uint32_t kSampleRate = 16000;  // Hz - the rate samples was captured at

FixedFFT<float, int16_t> fft;  // float calc type, int16_t input samples
int16_t samples[kN];           // e.g. a 16-bit ADC/I2S reading buffer

// Bin k covers the frequency k * kSampleRate / kN.
float binFrequency(uint16_t k) { return (float)k * kSampleRate / kN; }

void setup() {
  Serial.begin(115200);
  while (!Serial) {}

  if (!fft.begin(kN)) {
    Serial.println("fft.begin() failed - kN must be a power of two <= FixedFFT<float>::kMaxN");
    return;
  }

  for (uint16_t t = 0; t < kN; t++) {
    float v = 0.4f * cosf(2.0f * PI * 20 * t / kN);
    samples[t] = (int16_t)lroundf(v * 32767.0f);
  }

  // int16_t input rescaled into float's [-1, 1) convention.
  fft.loadReal(samples, kN);

  // For float, fft() always returns 0: its output is the plain
  // unscaled DFT sum (no per-stage 1/2 scaling needed - see the class
  // doc comment in FixedFFTEngine.h), so magnitude() values are
  // directly comparable to e.g. numpy's np.fft.fft() output.
  fft.fft();

  uint16_t peakBin = 0;
  float peakMag = 0;
  for (uint16_t k = 0; k < kN / 2; k++) {
    float m = fft.magnitude(k);
    if (m > peakMag) {
      peakMag = m;
      peakBin = k;
    }
  }
  Serial.print("peak bin: ");
  Serial.print(peakBin);
  Serial.print("\tfrequency(Hz): ");
  Serial.print(binFrequency(peakBin));
  Serial.print("\tmagnitude: ");
  Serial.println(peakMag);
}

void loop() {}
