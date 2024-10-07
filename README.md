# Foundamental limitations of Multiple-proposals Markov chain Monte Carlo algorithms
This repository is associated with the article [Pozza and Zanella (2024+). *On the fundamental limitations of multiple proposals Markov chain Monte Carlo algorithms*]. The **key contribution of this paper is outlined below**.

>Multiple-proposals schemes are a class of Markov chain Monte Carlo methods that generate, at each iteration, the next state of the chain by selecting it from several candidate points. In spite of their popularity, their theoretical properties are still not  well understood. We fill such a gap by first introducing a simple, yet powerful, stylized representation able to link every multiple-proposals scheme to a simpler single-proposal algorithm. This result is then exploited to study how the spectral gap of any multiple-proposals algorithm can increase with number of proposed points.

This repository provides **code to reproduce the simulation results presented in Section 4 of the paper**. To do so, create the folder `Logistic` and download the `R` codes below into it. Finally, run `Simulation logistic regression.R` after selecting `Logistic` as the working directory. See the corresponding R code for more information on the parameters of each function.

- [`Simulation logistic regression.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Simulation%20logistic%20regression.R). This file contains the code needed to reproduce the simulation study on high-dimensional logistic regression reported in Section 4 of the paper.
  
- [`MTM algorithm.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/MTM%20algorithm.R). R implementation for the multiple-try Metropolis algorithms with random walk proposal described in Section 4 of the paper.

- [`Generalised Metropolis.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Generalised%20Metropolis.R). R implementation for the generalised Metropolis algorithm with random walk proposal described in Section 4 of the paper.

- [`MTM MALA algorithm.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/MTM%20MALA%20algorithm.R). R implementation for the multiple-try Metropolis algorithms with MALA proposal described in Section 4 of the paper.

- [`MALA algorithm.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/MALA%20algorithm.R). R implementation for single-proposal MALA algorithm.

- [`Log_densities.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Log_densities.R). R code containing the log-densities for the posterior distribution of a Bayesian logistic model and for a multivariate normal distribution.

- [`Log_derivatives.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Log_derivatives.R). R code containing the first derivative for the log-density reported in `Log_densities.R`
- [`Summary stats.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Summary%20stats.R). R code to obtain the Expected Squared Jumping Distance (ESJD) of a Markov chain.

  

All the analyses are performed with a **MacBook Pro (OS Sonoma, version 14.5)**, using a `R` version **4.3.2**.

IMPORTANT: Although a seed is set at the beginning of the simulation for reproducibility, the final output may be subject to slight variations depending on which versions of the R packages were used to implement the code. This is due to possible internal changes to certain functions when the package version is updated. However, the magnitude of these small variations is negligible and does not affect the final conclusions.

