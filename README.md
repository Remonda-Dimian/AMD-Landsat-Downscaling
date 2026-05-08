# AMD-Landsat-Downscaling

This repository contains the key implementation scripts for the Adaptive Multi-Pixel Downscaling (AMD) framework developed for generating 10 m Landsat reflectance products across contrasting arid ecosystems.

## Repository Structure
* **GEE/**: Google Earth Engine (JavaScript) scripts for:
    * Landsat-7 SLC-off gap filling using Landsat-5
    * Adaptive Multi-Pixel Downscaling (AMD)
* **Validation/**: Python scripts for calculating:
    * Gradient Resolution Index (Grad. RI)
    * Frequency Resolution Index (Freq. RI)
    * Structural Similarity Index (SSIM. RI)
    * Weighted Resolution Index (Weighted RI)

## Prerequisites
To run the validation scripts, ensure you have Python 3.x installed along with the following libraries:
`numpy`, `rasterio`, `scipy`, `opencv-python`, `scikit-image`

## Sensors
* Landsat-7 ETM+
* Landsat-8 OLI
* Sentinel-2 MSI

## Study Areas
* Hotan, China (cold desert ecosystem)
* Kharga, Egypt (hot desert ecosystem)

## Temporal Coverage
* Landsat-7: 2000–2003, 2006–2010
* Landsat-8 & Sentinel-2: 2018–2024

## Method Summary
AMD integrates spectral similarity, spatial proximity, and texture-based adaptive weighting to reconstruct physically consistent 10 m reflectance imagery without external training datasets.

## Citation
If you use this repository, please cite the associated publication. [![DOI](https://zenodo.org/badge/1231517485.svg)](https://doi.org/10.5281/zenodo.20077458)

## License
This repository is distributed under the MIT License.
