# tensors

An early-stage personal repository for numerical-computing and physics
experiments. Its current contents are a single, self-contained project: a
real-time 2D wave interference and diffraction simulator for macOS.

## Status

This repository is a work in progress. At present the only implemented project
is the slit-diffraction simulator described below, under `wave-cycles/`. There
is no Rust or machine-learning code here yet, despite the repository name.

## Contents

- `wave-cycles/` — a real-time 2D wave interference and diffraction simulator
  (macOS, Cocoa, Objective-C++).

## Slit-diffraction simulator

An interactive Cocoa application that visualizes single- and multi-slit
diffraction in real time. The wavefield is computed directly from the
Huygens-Fresnel principle: each slit is discretized into point sources, and the
field at every pixel is the superposition of 2D cylindrical waves from those
sources. The visible wavelength is mapped to an approximate RGB color, and the
field computation is parallelized across pixel rows using Grand Central Dispatch.

Interactive controls:

- **Wavelength** (380-780 nm)
- **Slit width**
- **Slit separation**
- **Number of slits** (1-5; single-slit diffraction or multi-slit interference)
- **Wave speed** with play/pause for time animation

### Build and run

Building natively requires macOS and the Cocoa framework. From `wave-cycles/`:

```bash
make
./SlitDiffraction
```

To compile off macOS you would need a cross-platform Objective-C runtime that
provides the Cocoa API (such as GNUstep).

## Authors

The diffraction simulator was written by Umberto Puddu and Risa Charvi Metta.
