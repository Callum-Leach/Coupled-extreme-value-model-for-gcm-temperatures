# Coupled Extreme Value Model For GCM Temperatures

## Downloaded Data
- [X] UKESM1-0-LL - [1, 2, 3, 4, 8] (note this will include f2 as there is not option for f1)
- [X] ACCESS-CM2 - [1, 2, 3, 4, 5]
- [X] MRI-ESM2-0 - [1, 2, 3, 4, 5], missing data was downloaded manually. 
- [X] CESM2 - [4, 10, 11]
- [X] EC-Earth3 - [1, 4, 11] (The r9 ensemble did not want to be downloaded so we have removed this ensemble in its entirety. Also remove r4f2)

## Locations
- AN (Antarctic) [-180, -90, 180, -60]
- DA (Lut Desert) [57, 28, 59.5, 32]
- MO (Mojave Desert) [-118.5, 34, -114, 37.25]
- SA (Sahharah Desert) [-17.0722, 9.341, 38.4614, 34.6087] (only 832 processed?)
- SI (Simpsons Desert) [133, -26, 138, -16]
- UK (UK) [-13.6913, 49.9096, 1.7712,, 60.8479]

## Appendix / Supplementary Material (SM)

- **Figure SM1: 6 x 6 time-series plot of means**
  - Follows the same structure as Figure 2
- **Figure SM2: 6 x 6 time-series plot of medians**
  - Follows the same structure as Figure 2
- **Figure SM3: 2 x 6 linear regression slopes**
  - Row 1: Mean
  - Row 2: Median

## Data Processing

As previous work except ...
Put in details of data processing from phil's email. XXX

## Story of paper
- model selection.
- Non-linearity of model parameters. 