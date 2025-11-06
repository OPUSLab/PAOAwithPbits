# PAOA with p-Bits
This repository is to generate and plot the results of the paper: Generalized Probabilistic Approximate Optimization Algorithm

The paper is available on arXiv here: https://arxiv.org/abs/2507.07420v2

This repository contains code and data to reproduce the results of the paper titled:  
**“Generalized Probabilistic Approximate Optimization Algorithm.”**

Each folder corresponds to a main figure in the paper, including code for training, inference, and plotting.

## ⚙️ Figure 2: Majority Gate
 
In this section, we use a **private inverse temperature schedule** (β) to learn the weights of a majority gate.

- Numerical training code is provided to reproduce the results.  
- Plotting scripts are also included to visualize the learned behavior.


## ⚙️ Figure 3: Simulated Annealing Discovery

Here, we apply our **online annealing architecture** to **discover Simulated Annealing** within the PAOA framework using **single-parameter schedule optimization**.

- Training, inference, and plotting scripts are available.  
- The setup demonstrates that PAOA can naturally discover simulated annealing dynamics.


## ⚙️ Figure 4: SK Model

This experiment applies PAOA with **dual annealing schedule parametrization** to the **Sherrington-Kirkpatrick (SK) spin glass model**.
- The PAOA is applied to 26-spin SK model using matched count parameters as QAOA.
- We run QAOA using optimized parameters from Farhi et al. for comparison.  
- Includes training, inference, problem instances, and visualization code.

## ⚙️ Figure 5: SK Model with Lévy Bonds

We extend the dual-schedule PAOA to the SK model with **Lévy-distributed couplings**, introducing **schedule heterogeneity based on bond strength**.

- Heavily connected nodes are assigned a **lower annealing profile**.  
- Weakly connected nodes are assigned a **higher annealing profile**.  
- Full training, inference, and plotting tools are provided.

## Contributing

Contributions to improve the code or extend its functionality are welcome. Please feel free to submit issues or pull requests.


## Acknowledgements

ASA, SC, and KYC acknowledge support from the
National Science Foundation (NSF) under award number
2311295, and the Office of Naval Research (ONR),
Multidisciplinary University Research Initiative (MURI)
under Grant No. N000142312708. We are grateful to
Navid Anjum Aadit for discussions related to the hardware
implementation of online annealing, and Ruslan Shaydulin
and Zichang He for input on QAOA benchmarking. Use
was made of computational facilities purchased with funds
from the National Science Foundation (CNS-1725797) and
administered by the Center for Scientific Computing (CSC). The CSC is supported by the California NanoSystems
Institute and the Materials Research Science and Engineering
Center (MRSEC; NSF DMR 2308708) at UC Santa Barbara


## Contact

If you have any questions or suggestions, please open an issue in this repository or contact Abdelrahman Abdelrahman (abdelrahman@ucsb.edu).


