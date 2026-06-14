# GGMplus-processing-SouthernThailand
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20687627.svg)](https://doi.org/10.5281/zenodo.20687627)

MATLAB processing pipeline for satellite-derived gravity data analysis applied to geothermal exploration in southern Thailand.

---

## Overview 
This repository contains MATLAB scripts and processed data used to investigate geothermal potential in Thailand through satellite-derived gravity analysis. The primary data source is the Global Gravity Model plus 2013 (GGMplus).

---

## Repository Structure
```
GGMplus-processing-SouthernThailand/
│
├── derivative.m              # Gravity derivatives & Total Horizontal Derivative (THD)
├── preprocessing.m           # Complete Bouguer anomaly from free-air anomaly
├── terraincorrection.m       # Terrain correction (inner zone wedge + outer zone prism/ring)
├── upwardcontinuation.m      # Upward continuation for regional/residual separation
│
└── data/
    ├── *.tif                 # GeoTIFF raster outputs (gravity anomaly maps)
    └── *.mat                 # MATLAB workspace files (intermediate & final products)
```

---

## Requirements
- ***MATLAB*** R2023a
- [inpaint_nans](https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans)-some interpolation steps use the inpaint_nans function by John D'Errico (D'Errico, 2025)
