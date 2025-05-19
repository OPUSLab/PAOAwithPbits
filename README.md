# PAOA with p-Computer

This repository contains code and data to reproduce the results of the paper titled:  
**“Probabilistic Approximate Optimization Algorithm (PAOA) with Probabilistic Computers.”**

Each folder corresponds to a main figure in the paper, including code for training, inference, and plotting.

## 🟦 Figure 2: Majority Gate

In this section, we use a **private inverse temperature schedule** (β) to learn the weights of a majority gate.

- Numerical training code is provided to reproduce the results.  
- Plotting scripts are also included to visualize the learned behavior.


## 🟦 Figure 3: Simulated Annealing Discovery

Here, we apply our **online annealing architecture** to **discover Simulated Annealing** within the PAOA framework using **single-parameter schedule optimization**.

- Training, inference, and plotting scripts are available.  
- The setup demonstrates that PAOA can naturally discover simulated annealing dynamics.


## 🟦 Figure 4: SK Model

This experiment applies PAOA with **dual annealing schedule parametrization** to the **Sherrington-Kirkpatrick (SK) spin glass model**.

- We also run QAOA using optimized parameters from Farhi et al. for comparison.  
- Includes training, inference, problem instances, and visualization code.
- 

## 🟦 Figure 5: SK Model with Lévy Bonds

We extend the dual-schedule PAOA to the SK model with **Lévy-distributed couplings**, introducing **schedule heterogeneity based on bond strength**.

- Heavily connected nodes are assigned a **lower annealing profile**.  
- Weakly connected nodes are assigned a **higher annealing profile**.  
- Full training, inference, and plotting tools are provided.
