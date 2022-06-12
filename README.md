# Technical comment on "From white to green: Snow cover loss and increased vegetation productivity in the European Alps"

This repository includes the scripts used to compute the analysis for the technical comment adressed to Rumpf et al. (2022) published in Science : https://www.science.org/doi/10.1126/science.abn6697

This repository includes two files :
- The Google Earth Engine script based on the script of Rumpf et al. (2022) but modified to extract the 'Image Count' more efficiently.
- The R code used to conduct the analysis and produce the figures.

----- Methods :

We downloaded the time series of 'Vegetation productivity', 'Summer snow', and 'Permanent snow' from 1984 to 2021 here : https://zenodo.org/record/6386268
We reworked the Google Earth Engine script from Rumpf et al. (2022) to simplify the extraction of 'Image count' data from 1984 to 2021. The new code can be found in this repository.
We reprojected all data to the ETRS89 (EPSG=3035) coordinate system and resampled pixels to 300 m (factor 13.3) to reduce the computation cost and duration.
We used the EU-DEM v1.1 (https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1) and reprojected it on the 300 m grid.


