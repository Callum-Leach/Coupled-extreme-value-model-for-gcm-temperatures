import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import genextreme
from scipy.linalg import sqrtm

from scipy.stats import genextreme as gev

class GevMCMC:

    def __init__(self, data, param_setup, verbose=True):

        """
        This function sets up the data and the model parameters and hyperparameters for fitting a GEV model.

        Input
        data: T x 3 DataFrame where each column denotes each scenario 
        """
        
        # self.data = data[~np.isnan(data), :]
        self.data = data.dropna()
        self.param_setup = param_setup
        self.verbose = verbose
        self.nT = len(self.data.iloc[:, 0])
        self.Tim = np.linspace(0, 1, self.nT)
        self.nPrm, self.param_structure = self.setup_params()

        self.initial_params = self.find_starting_parameters()

    def setup_params(self):
        """
        This function sets up the correct parameters given a string. It sets up the total number of parameters, along with the parameter indices.

        Parameter arrangement is
         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  (21)
         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16      17      18      (18)
         1   2   3   4   5   6   7   8   9      10      11      12  13      14      15      (15)
         1   2       3       4       5   6       7       8       9  10      11      12      (12)
         1   2       3       4       5   6       7       8       9                          ( 9)
         1   2       3       4       5                           6                          ( 6)
         1                           2                           3                          ( 3)
        m0 m11 m12 m21 m22 m31 m32  s0 s11 s12 s21 s22 s31 s32  x0 x11 x12 x21 x22 x31 x32
        """

        config_map = {
        "CCC": (3, [0, 0, 0]),       # Only constant terms
        "LCC": (6, [1, 0, 0]),       # Linear Mu, constant Sgm and Xi
        "LLC": (9, [1, 1, 0]),       # Linear Mu and Sgm, constant Xi
        "QCC": (9, [2, 0, 0]),       # Quadratic Mu, constant Sgm and Xi
        "LLL": (12, [1, 1, 1]),      # Linear for all
        "QLC": (12, [2, 1, 0]),      # Quadratic Mu, Linear Sgm, constant Xi
        "QLL": (15, [2, 1, 1]),      # Quadratic Mu, Linear Sgm and Xi
        "QQC": (15, [2, 2, 0]),      # Quadratic Mu and Sgm, constant Xi
        "QQL": (18, [2, 2, 1]),      # Quadratic Mu and Sgm, Linear Xi
        "QQQ": (21, [2, 2, 2]),      # Quadratic MU, Sgm, Xi
        }

        if self.param_setup not in config_map:
            raise ValueError(f"Invalid parameter setup: {self.param_setup}")
        
        nPrm, param_structure = config_map[self.param_setup]

        return nPrm, param_structure
    
    def setup_dataframes(self):
        """
        Dynamically sets up the `samples` and `total_accepted` DataFrames 
        based on the parameter structure and number of scenarios.
        """
        column_names = []
        for param_type, param_name in zip(self.param_structure, ["mu", "s", "x"]):
            if param_type >= 0:  # Constant term
                column_names.append(f"{param_name}_0")
            for degree in range(1, param_type + 1):  # Loop through degrees (1 = Linear, 2 = Quadratic)
                for scenario in range(1, 4):  # 3 scenarios
                    column_names.append(f"{param_name}_{degree}{scenario}")

        # Create the DataFrames with the generated columns
        samples = pd.DataFrame(columns=column_names)
        total_accepted = pd.DataFrame(columns=column_names)

        return samples, total_accepted
    
    def build_parameter_arrays(self, params, time_steps):

        Mu = np.zeros((self.nT, 3))
        Sgm = np.zeros((self.nT, 3))
        Xi = np.zeros((self.nT, 3))

        # Compute Mu, Sgm, and Xi for each scenario
        for iS in range(3):

            """
            # Dynamic starting positions:
            Mu = 0
            Sgm = Mu.degree*3 + 1
            Xi = 3*(Mu.degree + Sgm.degree) + 2

            So what is the best way forward. We firstly want to iterate over all scenarios.
            Then look at the degree of each parameter
            """

            Mu_start = 0
            Sgm_start = 3 * self.param_structure[0] + 1
            Xi_start = 3 * (self.param_structure[0] + self.param_structure[1]) + 2

            # Calculate Mu for this scenario
            Mu[:, iS] = params[Mu_start]  # Constant term
            if self.param_structure[0] >= 1:  # Linear term
                Mu[:, iS] += params[Mu_start + iS + 1] * time_steps
            if self.param_structure[0] >= 2:  # Quadratic term
                Mu[:, iS] += params[Mu_start + 4 + iS] * (time_steps ** 2)

            # Calculate Sgm for this scenario
            Sgm[:, iS] = params[Sgm_start]  # Constant term
            if self.param_structure[1] >= 1:  # Linear term
                Sgm[:, iS] += params[Sgm_start + iS + 1] * time_steps
            if self.param_structure[1] >= 2:  # Quadratic term
                Sgm[:, iS] += params[Sgm_start + 4 + iS] * (time_steps ** 2)

            # Calculate Xi for this scenario
            Xi[:, iS] = params[Xi_start]  # Constant term
            if self.param_structure[2] >= 1:  # Linear term
                Xi[:, iS] += params[Xi_start + iS + 1] * time_steps
            if self.param_structure[2] >= 2:  # Quadratic term
                Xi[:, iS] += params[Xi_start + 4 + iS] * (time_steps ** 2)

        return Mu, Sgm, Xi
        
    def log_likelihood(self, params, time_steps):
        """
        This function calculates the log_likelihood of a generalised extreme value distribution

        Input (list, list): params which is a list of estimates for the three gev parameters. time_steps is an array of points from 0,1 which is the length of the data.

        Output (float): the log_likelihood of a generalised extreme value distribution with the passed in parameters.
        """

        tNll = np.zeros([3, 1])

        Mu, Sgm, Xi = self.build_parameter_arrays(params, time_steps)

        for iS in range(3):

            x = self.data.iloc[:, iS]

            # Compute the standardized variable (z)
            z = (x - Mu[:, iS]) / Sgm[:, iS]
            if np.any(1 + Xi[:, iS] * z <= 0):
                return np.inf  # Invalid parameters (domain constraint violated)

            # Compute the log-likelihood for GEV
            if np.abs(Xi[:, iS]).max() < 1e-6:  # Limiting case as xi -> 0 (Gumbel)
                log_pdf = -np.log(Sgm[:, iS]) - z - np.exp(-z)
            else:  # General case for xi != 0
                t = 1 + Xi[:, iS] * z
                log_pdf = (
                    -np.log(Sgm[:, iS])
                    - (1 / Xi[:, iS] + 1) * np.log1p(Xi[:, iS] * z)
                    - t**(-1 / Xi[:, iS])
                )

            # Accumulate the negative log-likelihood
            tNll[iS] = -1 * np.sum(log_pdf)

        return np.sum(tNll)

    
    def log_prior(self, params, time_steps):
        """
        This function provides prior constraints on estimated parameters.

        Input (list, list): params which is a list of estimates for the three gev parameters. time_steps is an array of points from 0,1 which is the length of the data.

        Output (float): Will return either 0 or -inf depending on the input parameters.
        """

        # m_0, m_11, m_21, m_31, s_0, s_11, s_21, s_31, x_0, x_11, x_21, x_31 = params

        Mu, Sgm, Xi = self.build_parameter_arrays(params, time_steps)

        for iS in range(3):
            
            t0=(1+Xi[:, iS]*(self.data.iloc[:, iS] - Mu[:, iS])/Sgm[:, iS])

            if (min(Sgm[:, iS]) <= 0) or (min(Xi[:, iS]) <= -1) or (max(Xi[:, iS]) > 0.2) or (min(t0) < 0):
                return -np.inf
        return 0
            
    
    def propose(
        self, params, iteration, beta, total_accepted, samples, burn_in, time_steps
    ):
        # Make sure the total_accepted array has float values
        # accepted_chain = np.array(total_accepted.iloc[max(0, iteration-999):, :], dtype=float)
        accepted_chain = np.array(total_accepted.iloc[max(0, len(total_accepted)-1000):, :], dtype=float)


        if iteration <= burn_in:
            new_parameters = [param + np.random.normal(0, 1) * 0.1 for param in params]

        else:
            SH = sqrtm(np.cov(accepted_chain, rowvar=False))
            SH = np.real(np.array(SH))
            z1 = np.random.normal(0, 1, size=self.nPrm)
            z2 = np.random.normal(0, 1, size=self.nPrm)
            y1 = (2.38 / np.sqrt(self.nPrm)) * np.matmul(SH, z1)
            y2 = (0.1 / np.sqrt(self.nPrm)) * z2
            new_parameters = params + (1 - beta) * y1 + beta * y2

        return new_parameters
    
    def acceptance_prob(self, old_params, new_params, time_steps):
        log_prior_old = self.log_prior(old_params, time_steps)
        log_prior_new = self.log_prior(new_params, time_steps)
        log_likelihood_old = self.log_likelihood(old_params, time_steps)
        log_likelihood_new = self.log_likelihood(new_params, time_steps)

        log_ratio = (-1*log_likelihood_new + log_prior_new) - (
            -1*log_likelihood_old + log_prior_old
        )

        return np.exp(log_ratio), log_likelihood_old

    def metropolis_hastings(
        self, n_samples, n2plt, burn_in=1000, thinning=10, beta=0.5, NGTSTR=0.1
    ):
        
        samples, total_accepted = self.setup_dataframes()
        ar = pd.DataFrame(columns=["acceptance_rate"])
        nloglikelihood = pd.DataFrame(columns=["negative_log_likelihood"])
        current_params = self.initial_params
        time_steps = np.linspace(0, 1, len(self.data.iloc[:, 0]))

        for i in range(n_samples):
            if i <= burn_in:
                for j in range(self.nPrm):
                    new_params = current_params.copy()

                    new_params[j] = new_params[j] + np.random.normal(0, 1) * NGTSTR

                    acceptance_probability, nll = self.acceptance_prob(
                        current_params, new_params, time_steps
                    )

                    if np.random.uniform(0, 1) < acceptance_probability:

                        current_params = new_params.copy()

                nloglikelihood = pd.concat(
                    [nloglikelihood, pd.DataFrame({"negative_log_likelihood": [nll]})],
                    ignore_index=True,
                )

                total_accepted = pd.concat(
                    [
                        total_accepted,
                        pd.DataFrame(
                            [current_params], columns=total_accepted.columns
                        ),
                    ],
                    ignore_index=True,
                )

            else:
                
                new_params = self.propose(
                    current_params,
                    i,
                    beta,
                    total_accepted,
                    samples,
                    burn_in,
                    time_steps,
                )
                
                acceptance_probability, nll = self.acceptance_prob(
                    current_params, new_params, time_steps
                )

                if np.random.uniform(0, 1) < acceptance_probability:

                    current_params = new_params.copy()
                
                nloglikelihood = pd.concat(
                    [nloglikelihood, pd.DataFrame({"negative_log_likelihood": [nll]})],
                    ignore_index=True,
                )

                total_accepted = pd.concat(
                    [
                        total_accepted,
                        pd.DataFrame(
                            [current_params], columns=total_accepted.columns
                        ),
                    ],
                    ignore_index=True,
                )

            if i >= n2plt and i % thinning == 0:
                samples = pd.concat(
                    [samples, pd.DataFrame([current_params], columns=samples.columns)],
                    ignore_index=True,
                )

            acceptance_rate = len(total_accepted) / (i + 1)

            ar = pd.concat(
                [ar, pd.DataFrame([acceptance_rate], columns=["acceptance_rate"])],
                ignore_index=True,
            )

        return samples, total_accepted, ar, nloglikelihood

    def find_starting_parameters(self):

        """
        This function aims to find suitable starting parameters pre mcmc.

        INPUT:
        OUTPUT:
        """

        # Here we simply fit a gumble distribution, setting xi = 0
        combined_data = self.data.values.ravel()

        beta = (np.sqrt(6) * np.std(combined_data)) / np.pi
        mu = np.mean(combined_data)
        xi = 0

        PrmCns=np.zeros(self.nPrm)

        param_idx = 0

        # Loop over the parameter structure (constant, linear, quadratic terms)

        for param_type, default_value in zip(self.param_structure, [mu, beta, xi]):
            # Assign the constant term (degree 0)
            PrmCns[param_idx] = default_value
            param_idx += 1

            # Assign zeros for linear and quadratic terms across 3 scenarios
            for degree in range(1, param_type + 1):  # Iterate over linear and quadratic terms
                for _ in range(3):  # Three scenarios for each degree
                    PrmCns[param_idx] = 0  # Initialize higher-order terms to zero
                    param_idx += 1
        
        return PrmCns
    
    def plot_return_values(self, samples):
        RtrPrd = 100  # Return Period
        nRls = 250  # Number of samples to use.

        # Get the time_step value for 2025 & 2125
        t_start = (2025 - 2015) / (2100 - 2015)
        t_end = (2125 - 2015) / (2100 - 2015)

        tXi = samples.iloc[:, 4] + t_start * samples.iloc[:, 5]
        tSgm = samples.iloc[:, 2] + t_start * samples.iloc[:, 3]
        tMu = samples.iloc[:, 0] + t_start * samples.iloc[:, 1]

        # Return value calculations at start
        RV_Start = genextreme.ppf(1 - 1 / RtrPrd, c=tXi, loc=tMu, scale=tSgm)

        # Parameter estimates at end
        tXi_End = samples.iloc[:, 4] + t_end * samples.iloc[:, 5]
        tSgm_End = samples.iloc[:, 2] + t_end * samples.iloc[:, 3]
        tMu_End = samples.iloc[:, 0] + t_end * samples.iloc[:, 1]

        # Return value calculations at end
        RV_End = genextreme.ppf(1 - 1 / RtrPrd, c=tXi_End, loc=tMu_End, scale=tSgm_End)
        RV_Delta = RV_End - RV_Start
        # Summary statistics
        Prb = np.nanmean(RV_End > RV_Start)
        Cdf = np.column_stack((np.sort(RV_Start), np.sort(RV_End)))
        Qnt = np.percentile(RV_Start, [2.5, 50, 97.5]), np.percentile(
            RV_End, [2.5, 50, 97.5]
        )

        RV_Delta_Return = [
            RV_Delta,
            RV_Start,
            RV_End,
            Prb,
            np.nanpercentile(RV_Delta, 2.5),
            np.nanmedian(RV_Delta),
            np.nanpercentile(RV_Delta, 97.5),
        ]

        return RV_Delta_Return
    

    def run(self, n_samples, n2plt, burn_in=1000, thinning=1, beta=0.05, NGTSTR=0.1):
        """
        Note that the genextremefunction uses the convention for the sign of the shape
        given in the documentation https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.genextreme.html
        We compensate for this when displaying graphs, but keep it for computation purposes.
        """

        samples, total_accepted, ar, nloglikelihood = self.metropolis_hastings(
            n_samples, n2plt, burn_in, thinning, beta
        )

        return samples, total_accepted, ar, nloglikelihood


        




