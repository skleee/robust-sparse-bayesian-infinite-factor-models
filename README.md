# Robust Sparse Bayesian Infinite Factor Models
* This repository contains an Rcpp-based R program for massive covariance estimation, using [robust sparse Bayesian infinite factor model](https://arxiv.org/abs/2012.04315) (Lee and Lee, 2020).
* The model is an extension of [Bhattacharya & Dunson (2011)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3419391/).
* The code is built based on source code of R package [`infinitefactor`](https://rdrr.io/github/poworoznek/infinitefactor/).

## Contents
* `MGPmain_t.R` : main file for covariance estimation.
  * `MGPprior.R`, `mgp_funcs.cpp` : source code file for *sparse Bayesian infinite factor models* (Bhattacharya & Dunson, 2011)
  * `MGPprior_t.R`, `mgp_funcs_t.cpp` : source code file for *robust sparse Bayesian infinite factor models* (Lee and Lee, 2020)
## References
* Lee, J., & Lee, J. (2020). Robust Sparse Bayesian Infinite Factor Models. *arXiv preprint arXiv:2012.04315*.
* Bhattacharya, A., & Dunson, D. B. (2011). Sparse Bayesian infinite factor models. *Biometrika*, 291-306.
