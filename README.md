# Coupled Extreme Value Model For GCM Temperatures

## Todo

**20/02/2025**

<del>1.  Create a new repo and upload all the downloaded data for the 6 regions (Antarctic, Simpsons Desert, Sahara Desert, Dasht-e Lut Desert, Mojave Desert, UK (base case)). Save data for the annual min, max, mean and median with a view to only use the min and max data for now.</del>

1. Generate diagnostic plots for the data. 
   - For each climate combination (gcm, Tas) plot time series of the scenario (SSP126, green; SSP245, orange; SSP585, grey) along with ensembles representing the alpha level. After this could look at the convergence of the mcmc chain.
   - **Figures to produce:**
     - **Figure 1: Map**
       - global map similar to figure 1, containing all locations.
     - **Figure 2: 6*6 time-series plot of max**
       - *Rows:* 6 different regions
       - *Columns:* 6 different GCMs
       - Panel will have 3 colours (per scenario), and up to line styles (ensemble members)
     - **Figure 3: 6*6 time-series plot of min**
       - Following figure 2.
     - **Figure 4: 2*6 linear regression slopes**
       - *Row 1:* Max
       - *Row 2:* Min
       - *Columns:* 6 different regions
       - Panel will be 6 x-values (GCMs), and coloured dots (per scenario)
       - Analogous to F3 from previous paper
     - **Table 1:**
       - List the details of all data extracted

2. Create a toy example to check the DIC calculations from MCMC (Use Phil's MATLAB code as reference) are working correctly & are consistent.
   1. For heatmap: Have colours indept of column, then use Rankscore for colours. Swap around the x-axis. So blue on leading diagonal.
   2. Calculate Rankerror(Normalised Rankscore), Rangeerror(as is), Rvdifferror(Bias in RVDeltas).
   3. Calculate the RVDelta for 2125 and 2025 based on the true and suggested models. Generate heatmap similar to before, x-axis: true model, y-axis: fitted model look at bias in RVDelta (0 is perfect).
   4. Put in plots for the above XXX
   5. Put in description of DIC and the method used for model selection. XXX

3. As well as DIC use 76 - 10 points to train and test the models predictive performance.
   1. Take predictive likelihood from 250 mcmc sets of model parameters, then evaluate the likelihood of the data for all the sets of parameters for the last 10 values.
   2. Take the avg of these -ve likelihoods giving a measure of goodness of fit. 
   3. Compare this with DIC.

4. Run the DIC model selection MCMC on the real data and store the results as before.

5. (Phil) to look into slab-spike methods & Reversible jump.

Rankscore - The rank of the model compared to all others.
Rankerror - (RankOfModel - 1) / (RankOfWorstModel - 1)
Rangeerror - (DicOfModel - DicBestModel) / (DicWorstModel - DicBestModel)
Rvdifferror - Bias in RVDeltas for each model.

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

## MCMC Runs
- [X] MR
- [] UK
- [] AC
- [] CE
- [] EC

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