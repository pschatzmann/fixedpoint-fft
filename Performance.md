# Microcontroller FFT & IFFT Performance Benchmark (N=64)
*A comparative analysis of DSP execution performance across multiple MCU architectures.*

---

## 1. Executive Summary & Key Takeaways
This benchmark measures the time taken to compute an **$N=64$ Complex FFT and IFFT** across different data types (`int8_t`, `int16_t`, `int32_t`, `float`, and `double`). The processors span from vintage 8-bit AVRs to high-performance ARM Cortex-M7 cores.

### Core Architectural Insights:
* **The Single-Precision Float King:** **STM32H7** dominates overall performance due to its massive 480 MHz clock speed and dedicated double-precision hardware FPU. However, the **ESP32-S3** and **ESP32** punch far above their weight class in single-precision `float`, outperforming everything except the STM32H7.
* **ESP32-S3 vs. ESP32 Generation Leap:** The ESP32-S3 (Xtensa LX7) demonstrates a massive upgrade over the base ESP32 (Xtensa LX6). It reduces floating-point FFT latency from $75.14\,\mu\text{s}$ down to **$62.62\,\mu\text{s}$** (a ~17% speedup) and crushes integer IFFTs, with `int16_t` running in just **$107.96\,\mu\text{s}$** (~33% faster).
* **The Double-Precision Cliff:** Devices without a double-precision hardware FPU (e.g., ESP32/S3, STM32F411, UNO R4, RP2040) experience a massive performance drop (frequently 5x to 11x slower) when stepping from 32-bit `float` to 64-bit `double` because double-precision math must be emulated in software.
* **RP2040 Integer Power:** The **RP2040**'s integer processing is exceptionally fast—beating even the STM32F411 in 8-bit and 16-bit FFTs—but falls behind on floats due to the lack of an on-chip hardware FPU (relying instead on its optimized ROM floating-point libraries).

---

## 2. FFT Master Comparison Tables

### FFT Calculation Latency (lower is better)
Values are in microseconds ($\mu\text{s}$) per single $N=64$ FFT calculation.

| Microcontroller | Core Architecture | Clock Speed | int8_t | int16_t | int32_t | float (32-bit) | double (64-bit) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Arduino Nano** | AVR ATMega328P (8-bit) | 16 MHz | 4,144.62 | 9,798.18 | 29,160.16 | 17,165.54 | *N/A* |
| **UNO R4** | ARM Cortex-M4 (32-bit) | 48 MHz | 732.87 | 668.05 | 1,108.19 | 360.89 | 4,045.90 |
| **STM32F411** | ARM Cortex-M4 (32-bit) | 100 MHz | 318.39 | 284.30 | 356.52 | 151.90 | 1,737.21 |
| **RP2040** | ARM Cortex-M0+ (32-bit) | 133 MHz | 228.24 | 228.37 | 472.07 | 939.43 | 1,672.06 |
| **ESP32** | Tensilica Xtensa LX6 (32-bit)| 240 MHz | 255.05 | 250.49 | 398.47 | 75.14 | 756.14 |
| **ESP32-S3** | Tensilica Xtensa LX7 (32-bit)| 240 MHz | 199.44 | 190.14 | 361.17 | 62.62 | 702.66 |
| **STM32H7** | ARM Cortex-M7 (32-bit) | 480 MHz | **44.76** | **36.85** | **45.81** | **23.36** | **238.10** |

### FFT Throughput (higher is better)
Values represent total **FFT operations computed per second** (FFT/s).

| Microcontroller | int8_t | int16_t | int32_t | float | double |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Arduino Nano** | 241 | 102 | 34 | 58 | *N/A* |
| **UNO R4** | 1,364 | 1,497 | 902 | 2,771 | 247 |
| **STM32F411** | 3,141 | 3,517 | 2,805 | 6,583 | 576 |
| **RP2040** | 4,381 | 4,379 | 2,118 | 1,064 | 598 |
| **ESP32** | 3,921 | 3,992 | 2,510 | 13,308 | 1,323 |
| **ESP32-S3** | 5,014 | 5,259 | 2,769 | 15,971 | 1,423 |
| **STM32H7** | **22,341** | **27,141** | **21,827** | **42,808** | **4,200** |

---

## 3. IFFT Master Comparison Tables

### IFFT Calculation Latency (lower is better)
Values are in microseconds ($\mu\text{s}$) per single $N=64$ Inverse FFT calculation.

| Microcontroller | Core Architecture | Clock Speed | int8_t | int16_t | int32_t | float (32-bit) | double (64-bit) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Arduino Nano** | AVR ATMega328P (8-bit) | 16 MHz | 4,033.66 | 9,302.30 | 26,380.36 | 18,900.00 | *N/A* |
| **UNO R4** | ARM Cortex-M4 (32-bit) | 48 MHz | 741.59 | 652.91 | 1,084.73 | 384.70 | 4,415.17 |
| **STM32F411** | ARM Cortex-M4 (32-bit) | 100 MHz | 312.65 | 280.38 | 336.16 | 163.91 | 1,878.80 |
| **RP2040** | ARM Cortex-M0+ (32-bit) | 133 MHz | 219.27 | 224.29 | 464.54 | 967.96 | 1,808.77 |
| **ESP32** | Tensilica Xtensa LX6 (32-bit)| 240 MHz | 170.51 | 161.18 | 195.37 | 83.29 | 830.59 |
| **ESP32-S3** | Tensilica Xtensa LX7 (32-bit)| 240 MHz | 122.90 | 107.96 | 153.10 | 69.11 | 776.86 |
| **STM32H7** | ARM Cortex-M7 (32-bit) | 480 MHz | **43.65** | **36.79** | **39.06** | **24.32** | **258.70** |

### IFFT Throughput (higher is better)
Values represent total **IFFT operations computed per second** (IFFT/s), derived directly from calculation runtimes ($1,000,000 / \text{Latency}$).

| Microcontroller | int8_t | int16_t | int32_t | float | double |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Arduino Nano** | 247 | 107 | 37 | 52 | *N/A* |
| **UNO R4** | 1,348 | 1,531 | 921 | 2,599 | 226 |
| **STM32F411** | 3,198 | 3,566 | 2,974 | 6,100 | 532 |
| **RP2040** | 4,560 | 4,458 | 2,152 | 1,033 | 552 |
| **ESP32** | 5,864 | 6,204 | 5,118 | 12,006 | 1,203 |
| **ESP32-S3** | 8,137 | 9,263 | 6,532 | 14,470 | 1,287 |
| **STM32H7** | **22,910** | **27,181** | **25,601** | **41,118** | **3,865** |

---

## 4. Dedicated Device Analysis

### 🔴 Arduino Nano (AVR ATMega328P @ 16 MHz)
The legacy 8-bit architecture struggles immensely with DSP tasks.
* **Why it's slow:** It lacks a hardware multiplier for larger integers (32-bit multiplication must be emulated in software via dozens of assembly instructions) and has no FPU. 
* **Data-type anomaly:** Notice that `float` ($17,165\,\mu\text{s}$) is actually faster than `int32_t` ($29,160\,\mu\text{s}$). This is because the 32-bit float library uses 24-bit mantissa optimizations, whereas standard 32-bit integer division/multiplication on an 8-bit core requires broad software-based register shifting.

### 🔵 UNO R4 (Renesas RA4M1 Cortex-M4 @ 48 MHz)
A massive modernization for the basic Arduino form factor.
* **FPU Power:** The onboard Cortex-M4 single-precision FPU makes `float` math incredibly quick ($360.89\,\mu\text{s}$), running nearly **11x faster** than its `double` math ($4,045.90\,\mu\text{s}$), which must be calculated in software.

### 🟣 RP2040 (Cortex-M0+ @ 133 MHz)
The Raspberry Pi silicon lacks a hardware FPU but features incredibly fast integer execution.
* **Integer Strength:** It beats the STM32F411 at `int8_t` and `int16_t` processing because of its higher clock speed (133 MHz vs 100 MHz).
* **Smart Floating-Point ROM:** Lacking an FPU, the RP2040 relies on optimized floating-point routines burned directly into its boot ROM. This makes its software `float` ($939\,\mu\text{s}$) and `double` ($1,672\,\mu\text{s}$) performance respectable, though still significantly slower than hardware-accelerated chips.

### 🟢 STM32F411 (Cortex-M4 @ 100 MHz)
A balanced DSP powerhouse.
* Equipped with ARM's DSP instruction extensions and a single-precision FPU, it displays exceptionally uniform integer scaling and rapid `float` speeds ($151.90\,\mu\text{s}$).

### 🟠 ESP32 (Xtensa LX6 @ 240 MHz)
A consumer audio favorite.
* **Single-Precision Optimization:** Because of its high core clock speed (240 MHz) and optimized single-precision pipeline, its `float` FFT is blindingly fast at **$75.14\,\mu\text{s}$** (over 13,300 FFTs per second). It beats everything below the STM32H7 by a wide margin.
* **The IFFT Speedup:** Interestingly, the ESP32 calculates integer IFFTs significantly faster than regular FFTs (e.g., `int32_t` IFFT takes **$195.37\,\mu\text{s}$** compared to **$398.47\,\mu\text{s}$** for the FFT). This highlights high efficiency in compiler execution paths for the backward-scaling butterfly passes of the inverse algorithm.

### ⚡ ESP32-S3 (Xtensa LX7 @ 240 MHz)
The real measurements highlight the incredible processing benefits of the upgraded LX7 core.
* **Upgraded Instruction Set:** The S3 core features **PIE (Processor Instruction Extensions)** which provides hardware-level vector processing capability. 
* **The Float Champ:** At **$62.62\,\mu\text{s}$** for a floating-point FFT (nearly **16,000 computations per second**), the ESP32-S3 is a formidable chip for real-time edge audio processing, DSP filtering, and spectrum analysis.
* **Integer IFFT Beast:** The ESP32-S3 processes `int16_t` IFFT in just **$107.96\,\mu\text{s}$** (resulting in an outstanding **9,263 IFFTs per second**).

### 👑 STM32H7 (Cortex-M7 @ 480 MHz)
The absolute performance benchmark.
* **Hardware Double Precision:** The STM32H7 contains a dual-precision hardware FPU. It computes double-precision `double` values ($238.10\,\mu\text{s}$) faster than the Arduino Nano or RP2040 can calculate simple `int8_t` math!
* **Blazing Clock:** Running at 480 MHz with a superscalar L1 cache pipeline allows it to crank out **42,808 float FFTs every single second**.