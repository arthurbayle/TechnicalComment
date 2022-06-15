# Technical comment on "From white to green: Snow cover loss and increased vegetation productivity in the European Alps"

This repository includes the scripts used to compute the analysis for the technical comment adressed to Rumpf et al. (2022) published in Science : https://www.science.org/doi/10.1126/science.abn6697

This repository includes two files :
- The Google Earth Engine script based on the script of Rumpf et al. (2022) but modified to extract the 'Image Count' more efficiently.
- The R code used to conduct the analysis and produce the figures.


------------------------------------------------------------------- FIGURE 1, A,B,C,D,E -------------------------------------------------------------------

We downloaded the time series of 'Vegetation productivity', 'Summer snow', and 'Permanent snow' from 1984 to 2021 here : https://zenodo.org/record/6386268
We reworked the Google Earth Engine script from Rumpf et al. (2022) to simplify the extraction of 'Image count' data from 1984 to 2021. The new code can be found in this repository.
We reprojected all data to the ETRS89 (EPSG=3035) coordinate system and resampled pixels to 300 m (factor 13.3) to reduce the computation cost and duration.
We used the EU-DEM v1.1 (https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1) and reprojected it on the 300 m grid.

For Figure 1A, we divided the dataset in three periods of equal length to represent the change unequal numbers of images available throughout the time series AND on an elevational gradient, showing that the unequality worsen with increasing altitudes.
For Figure 1B, it is a figure similar to FIG S1b of Rumpf et al. (2022).
For Figure 1C, D, E, we randomly selected 100.000 pixels over the entire study area of Rumpf et al. (2022) for 'vegetation productivity', 'summer snow' and 'permanent snow'. We only showed results for images number between 1 and 10 for clarity as the relation is asymptotic.
For Figure 1D, we limit our analysis to pixels WITH snow, which we selected by removing pixels with a probability of summer snow inferior to 0.1. It was done as almost the entire dataset has a probability of 0 using the method of Rumpf et al. (2022).
For Figure 1E, as the variable 'permanent snow' is a binary variable, we computed the proportion of the 100.000 pixels that is classified as permenant snow and computed the relation with the number of images available per year.

------------------------------------------------------------------- FIGURE 1, F, G -------------------------------------------------------------------
