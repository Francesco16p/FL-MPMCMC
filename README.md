# Foundamental limitations of Multiple-proposals Markov chain Monte Carlo algorithms
This repository is associated with the article [Pozza and Zanella (2024+). *On the fundamental limitations of multiple proposals Markov chain Monte Carlo algorithms*]. The **key contribution of this paper is outlined below**.

>Multiple-proposals schemes are a class of Markov chain Monte Carlo methods that generate, at each iteration, the next state of the chain by selecting it from several candidate points. From a theoretical prospective, however, the effective gain provided by the adoption of these techniques is still not fully understood, mainly due to their apparent complexity and the differences which characterise alternative solutions within this class. This paper fill such a gap by first introducing a simple, yet powerful, stylized representation able to link every multiple-proposals scheme to a simpler single-proposal algorithm. Such a result is then exploited to study how the spectral gap of any multiple-proposals algorithm can increase with number of proposed points $K$. Such a theoretical investigation highlights that, although in general the improvement can depend at most linearly on $K$, in many practically relevant cases it is more limited and quantifiable in a poly-logarithmic gain.

This repository provides **code to implement the simulation results presented in the paper**. 

- [`CushingLogit.md`](https://github.com/Francesco16p/SMA/blob/main/CushingLogistic.md). 
  
- [`CushingProbit.md`](https://github.com/Francesco16p/SMA/blob/main/CushingLogistic.md). 

All the analyses are performed with a **MacBook Pro (OS Sonoma, version 14.5)**, using a `R` version **4.3.2**.

