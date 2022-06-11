# Technical comment on "From white to green: Snow cover loss and increased vegetation productivity in the European Alps"

This repository includes the scripts used to compute the analysis for the technical comment adressed to Rumpf et al. (2022) published in Science : https://www.science.org/doi/10.1126/science.abn6697

Arthur Bayle, Simon Gascoin, Bradley Z. Carlson, Christophe Randin, Edoardo Cremonese, Gianluca Filippa, Philippe Choler

This repository includes two files :
- The Google Earth Engine script based on the script of Rumpf et al. (2022) but modified to extract the 'Image Count' more efficiently.
- The R code used to conduct the analysis and produce the figures.

                                                          -------------- Figure 1. --------------

We downloaded the time series of 'Vegetation productivity', 'Summer snow', and 'Permanent snow' from 1984 to 2021 here : https://zenodo.org/record/6386268
We reworked the Google Earth Engine script from Rumpf et al. (2022) to simplify the extraction of 'Image count' data from 1984 to 2021. The new code can be found in this repository.
We reprojected all data to the ETRS89 (EPSG=3035) coordinate system and resampled pixels to 300 m (factor 10) to reduce the computation cost and duration.
We used the EU-DEM v1.1 (https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1) and reprojected it on the 300 m grid.

Figure 1B : We splitted the time series of 'Image count' for representation purposes and plotted the distribution of clear-sky pixels numbers according to elevation for the three period (1984-1996, 1997-2009 and 2010-2021)

Figure 1C : We plotted the distribution of clear-sky pixels numbers from 1984 to 2021 and showed the observation periods of the three satellites (L5, L7 and L8).

Figure 1D : We plotted the distribution of 'Vegetation productivity' (75% quantile of NDVI) for 1 to 10 clear-sky pixels numbers availability, showing a clear asymptotic relation that should be interpreted as the bias.
