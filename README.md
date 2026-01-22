# CMMRI: Color Sensor Accuracy Index

A MATLAB implementation of the **Camera Metamer Mismatch Radii Index (CMMRI)**. This metric evaluates the colorimetric accuracy of digital cameras by calculating the Metamer Mismatch Body (MMB) induced by the camera sensors relative to the human standard observer.

## Reference
```bibtex
@article{Roshan2020,
    author = {Roshan, Emitis and Funt, Brian},
    title = {Color Sensor Accuracy Index Utilizing Metamer Mismatch Radii},
    journal = {Sensors},
    volume = {20},
    number = {15},
    year = {2020},
    article-number = {4275},
    url = {https://www.mdpi.com/1424-8220/20/15/4275},
    doi = {10.3390/s20154275}
}

```

## Usage

The main entry point for the algorithm is `CMMRI.m`.

1. Ensure the data files `Canon500D.mat` and `D65_6500k.xlsx` are in your MATLAB path.
2. Run the script:

```matlab
CMMRI
```

The script calculates the Metamer Mismatch Volume (MMV) for a mid-grey sample and outputs the average MMV dimension.

## Key Files

* **`CMMRI.m`**: The main execution script.
* **`Canon500D.mat`**: Sample camera sensor sensitivities.
* **`D65_6500k.xlsx`**: Standard D65 illuminant data.
* **`calculate_mmv.m`**: Optimization routine for finding metameric boundaries.

## Requirements

* **MATLAB** (R2020a or later recommended)
* **Optimization Toolbox**
