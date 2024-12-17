# Coupled Extreme Value Model For GCM Temperatures

## Todo

**16/12/2024**

<del>1.  Create a new repo and upload all the downloaded data for the 6 regions (Antarctic, Simpsons Desert, Sahara Desert, Dasht-e Lut Desert, Mojave Desert, UK (base case)). Save data for the annual min, max, mean and median with a view to only use the min and max data for now.</del>

2. Generate diagnostic plots for the data. For each climate combination (gcm, Tas) plot time series of the scenario (SSP126, green; SSP245, orange; SSP585, grey) along with ensembles representing the alpha level. After this could look at the convergence of the mcmc chain.

3. Create a toy example to check the DIC calculations from MCMC (Use Phil's MATLAB code as reference) are working correctly & are consistent.

4. Run the DIC model selection mcmc on the real data and store the results as before.

5. (Phil) to look into slab-spike methods.

**09/12/2024**

Look at 6 climate regions:
- Sahara Desert
- Dasht-e Lut
- Mojave Desert
- Simpsons Desert
- Antarctic (Cold Stuff)
- UK (Baseline)

In order to choose the locations we can either do:

- Take a centre point, then choose closest 20 points to that centre.
- Choose a bounding box for each region, selecting all points within. 

Methodology:

- Continue with DIC and model selection.
- Use Reverse Jump MCMC