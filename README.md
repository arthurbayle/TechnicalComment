# Technical comment on "From white to green: Snow cover loss and increased vegetation productivity in the European Alps"

This repository includes the scripts used to compute the analysis for the technical comment adressed to Rumpf et al. (2022) published in Science : https://www.science.org/doi/10.1126/science.abn6697

This repository includes two files :
- The Google Earth Engine script based on the script of Rumpf et al. (2022) but modified to extract the 'Image Count' more efficiently.
- The R code used to conduct the analysis and produce the figures.


---------------------------------------------------------- FIGURE 1, A,B,C,D,E ---------------------------------------------------------

We downloaded the time series of 'Vegetation productivity', 'Summer snow', and 'Permanent snow' from 1984 to 2021 here : https://zenodo.org/record/6386268
We reworked the Google Earth Engine script from Rumpf et al. (2022) to simplify the extraction of 'Image count' data from 1984 to 2021. The new code can be found in this repository.
We reprojected all data to the ETRS89 (EPSG=3035) coordinate system and resampled pixels to 300 m (factor 13.3) to reduce the computation cost and duration.
We used the EU-DEM v1.1 (https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1) and reprojected it on the 300 m grid.

For Figure 1A, we divided the dataset in three periods of equal length to represent the increasing numbers of images available throughout the time series AND on an elevational gradient, showing that the inequality worsens with increasing elevation.
For Figure 1B, it is a figure similar to FIG S1b of Rumpf et al. (2022).
For Figure 1C, D, E, we randomly selected 100.000 pixels over the entire study area of Rumpf et al. (2022) for 'vegetation productivity', 'summer snow' and 'permanent snow'. We only showed results for images number between 1 and 10 for clarity as the relation is asymptotic.
For Figure 1D, we limit our analysis to pixels WITH snow, which we selected by removing pixels with a probability of summer snow inferior to 0.1. It was done as almost the entire dataset has a probability of 0 using the method of Rumpf et al. (2022).
For Figure 1E, as the variable 'permanent snow' is a binary variable, we computed the proportion of the 100.000 pixels that is classified as permanent snow and computed the relation with the number of images available per year.

![FIGURE1_Part1](https://user-images.githubusercontent.com/77186981/173824375-b5cf6dc2-cbf4-44c3-bed4-36a30bf39737.png)

------------------------------------------------------------ FIGURE 1, F, G ------------------------------------------------------

To estimate the bias highlighted in this technical comment, we downloaded each image of the time series for the Path/Row 195028 only. Our example is thus quantitatively representative of areas where there is no overlapping of two tiles which increases the number of images available. In case of overlapping tiles, the bias is expected to be lower as the bias follows an asymptotic relation.

The analysis was done on an area ranging from the Mont-Blanc massif to Zermatt as shown below :
![image](https://user-images.githubusercontent.com/77186981/173767392-4b08bbcc-7b11-4370-ab05-a98bd3fca591.png)


We randomly selected 10.000 pixels from pixels that are vegetated (NDVI 0.75 quantile superior to 0 on the image from 2021-07-29) over this area.

Greening trends are computed using the Theil-sen estimator (package 'mblm') and significance using the Mann-Kendall test (package 'Kendall'). 'DescTools' and 'zoo' packages were also used. Trends were computed only when 12 year of relevant data was available as done in Rumpf et al. (2021).

We first computed the greening trends using all available images (hereafter referred as ALL) as a basis.
Then we computed the greening trends 50 times by randomly selecting a maximum of 3 images per year (hereafter referred as 3MAX). For example, there are 15 images available in 2021 (8 from Landsat 7 and 7 from Landsat 8). For each iteration, we randomly 3 images only among the 15 available. In 1984, only 2 images were available. So for each iteration, the two same images were selected. So only years with more than 3 images can be considered as randomized.

When considering the average number of clear-sky pixels available per year, our randomization has resulted in 1.8 clear-sky pixels used per year.
![FIGURE1F](https://user-images.githubusercontent.com/77186981/173769170-2b95230e-a5a3-47a5-a2f2-7ed4ba5ae91b.png)

We calculated the proportion of pixels that experienced significant greening (p-val <0.05 and <0.01) based on the ALL and 3MAX methods. We found large difference in the % of pixels experiencing greening. This analysis is purely preliminary with the sole intention to portray the bias rather than trully estimates the "true" proportion of greening pixels. Indeed, using our method, we did not correct for the bias but simply remove it by degrading our data. See Berner et al. (2020) for methods to correct for the bias.

![FIGURE1G](https://user-images.githubusercontent.com/77186981/173769846-c4fd59d5-46ec-42c6-8f47-14321fea6842.png)
