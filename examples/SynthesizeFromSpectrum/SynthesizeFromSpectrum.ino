// Synthesizes a real time-domain signal from a frequency-domain
// specification, using loadHalfSpectrum() + ifft().
//
// ifft() only produces a purely real result (imaginary part ~0 at
// every sample) when the spectrum has Hermitian symmetry: X[N-k] =
// conj(X[k]). loadHalfSpectrum() takes just the non-redundant half of
// the spectrum - bins 0..N/2 - and mirrors it into the upper half
// automatically, so you only have to specify N/2+1 complex values
// instead of building the full N-length symmetric spectrum by hand.
//
// This sketch builds a spectrum with two frequency components set
// directly (a simple additive synthesizer) and inverse-transforms it
// into a time-domain waveform.

#include "FixedFFT.h"

using namespace fixedpoint_fft;

constexpr uint16_t kN = 64;
constexpr uint16_t kHalfN = kN / 2 + 1;

FixedFFT<int16_t, int16_t> fft;  // int16_t (Q15) calc type, int16_t input samples

void setup() {
  Serial.begin(115200);
  while (!Serial) {}

  if (!fft.begin(kN)) {
    Serial.println("fft.begin() failed");
    return;
  }

  // Build a one-sided spectrum: two bins with nonzero magnitude, real
  // part only (imag left at 0 for both - a pure cosine at each
  // frequency, no phase shift).
  int16_t real_in[kHalfN] = {0};
  int16_t imag_in[kHalfN] = {0};
  real_in[4] = 12000;   // stronger component at bin 4
  real_in[10] = 6000;   // weaker component at bin 10

  fft.loadHalfSpectrum(real_in, imag_in, kHalfN);

  // This spectrum didn't come from fft(), so it has no guaranteed
  // headroom - use the adaptive overflow guard (see ifft()'s doc
  // comment in FixedFFTEngine.h).
  fft.ifft(/*safeScale=*/true);

  Serial.println("sample\tvalue");
  for (uint16_t t = 0; t < kN; t++) {
    Serial.print(t);
    Serial.print('\t');
    Serial.println(fft.real(t));  // imag(t) is ~0 at every sample
  }
}

void loop() {}
