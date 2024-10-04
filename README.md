# Foundamental limitations of Multiple-proposals Markov chain Monte Carlo algorithms
This repository is associated with the article [Pozza and Zanella (2024+). *On the fundamental limitations of multiple proposals Markov chain Monte Carlo algorithms*]. The **key contribution of this paper is outlined below**.

>Multiple-proposals schemes are a class of Markov chain Monte Carlo methods that generate, at each iteration, the next state of the chain by selecting it from several candidate points. From a theoretical prospective, however, the effective gain provided by the adoption of these techniques is still not fully understood, mainly due to their apparent complexity and the differences which characterise alternative solutions within this class. We fill such a gap by first introducing a simple, yet powerful, stylized representation able to link every multiple-proposals scheme to a simpler single-proposal algorithm. This result is then exploited to study how the spectral gap of any multiple-proposals algorithm can increase with number of proposed points $K$. Such a theoretical investigation highlights that, although in general the improvement can depend at most linearly on $K$, in many practically relevant cases it is more limited and quantifiable in a poly-logarithmic gain.

This repository provides **code to implement the simulation results presented in Section 4 the paper**. To reproduce the results of the simulation study, create the folder `Logistic` and download the `R` codes below into it. Finally, run `Simulation logistic regression.R` after selecting the `Logistic` as the working directory. For more information on the parameters of each of the functions used in the sumulation study, see the corresponding R code. 

- [`Simulation logistic regression.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Simulation%20logistic%20regression.R). This file contains the code needed to reproduce the simulation study on high-dimensional logistic regression reported in Section 4 of the paper.
  
- [`MTM algorithm.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/MTM%20algorithm.R). R implementation for the multiple-try Metropolis algorithms with random walk proposal described in Section 4 of the paper.

- [`Generalised Metropolis.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Generalised%20Metropolis.R). R implementation for the generalised Metropolis algorithm with random walk proposal described in Section 4 of the paper.

- [`MTM MALA algorithm.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/MTM%20MALA%20algorithm.R). R implementation for the multiple-try Metropolis algorithms with MALA proposal described in Section 4 of the paper.

- [`MALA algorithm.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/MALA%20algorithm.R). R implementation for single-proposal MALA algorithm.

- [`Log_densities.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Log_densities.R). R code containing the log-densities for the posterior distribution of a Bayesian logistic model and for a multivariate normal distribution.

- [`Log_derivatives.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Log_derivatives.R). R code containing the first derivative for the log-density reported in `Log_densities.R`
- [`Summary stats.R`](https://github.com/Francesco16p/FL-MPMCMC/blob/main/Summary%20stats.R). R code to obtain Expected Squared Jumping Distance (ESJD) of a Markov chain.

  

All the analyses are performed with a **MacBook Pro (OS Sonoma, version 14.5)**, using a `R` version **4.3.2**.

IMPORTANT: Although a seed is set at the beginning of the simulation for reproducibility, the final output may be subject to slight variations depending on which versions of the R packages were used to implement the code. This is due to possible internal changes to certain functions when the package version is updated. However, the magnitude of these small variations is negligible and does not affect the final conclusions.

