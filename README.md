# ECG-QRS-detector
MATLAB implementation of the Pan-Tompkins QRS detection algorithm with adaptive thresholding (SPKI/NPKI), notch filtering, and spectrogram analysis for ECG signal processing.

> Digital Signal Processing project completed during a 15-hour course at 
> [ENSEIRB-MATMECA](https://enseirb-matmeca.bordeaux-inp.fr/), Bordeaux INP — 
> Telecom department, 2025/2026.

## Features
- Band-pass filtering (5–15 Hz)
- Notch filter at 12.5 Hz (power line noise)
- Adaptive SPKI/NPKI threshold (baseline drift robustness)
- Spectrogram analysis
- Fibrillation detection

## Files
| File | Description |
|------|-------------|
| `detect_R.m` | Main R-peak detection pipeline (PQRST + heart frequency) |
| `spectro_1.m` | Short-time Fourier transform spectrogram |
| `fibrillation.m` | Pathological ECG visualization |

## Usage
1. Load a `.mat` ECG file (variables: `ecg`, `Fs`)
2. Run `detect_R.m`

## Reference
Pan & Tompkins, *A real-time QRS detection algorithm*,
IEEE Trans. Biomed. Eng., 1985.
