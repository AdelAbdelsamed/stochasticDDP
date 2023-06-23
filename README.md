# Stochastic DDP

This project was conducted within the Seminar Learning, Control and Identification Algorithms during the summer semester of 2023 at the Chair of Intelligent Control Systems, RWTH (Rheinisch-Westfälische Technische Hochschule).
The seminar focused on *Stochastic Differential Dynamic Programming* (SDDP) based on the paper by Theodorou et al. [1]. The proposed algorithm was implemeneted in MATLAB. Implementational issues such as line-search and regularization are explored. The algorithm is demonstrated on three models with varying degrees of complexity,
including a 7-DOF parafoil homing trajectory optimization problem. 


## Table of Contents

- [Project Overview](#project-overview)
- [Usage](#usage)
- [Acknowledgement](#acknowledgement)
- [References](#references)


## Project Overview

The SDDP function can be applied to all Optimal Control Problems and must not be changed. Eventually, one can tune the hyperparameters of the SDDP algorithm. However, what changes with the OCP at hand are the cost function and the dynamics. 
Therefore, these are defined as abstract classes in `Cost_Fcn.mat` and `Dynamics.mat`. To define the OCP, two classes must be defined: one that inherits `Cost_Fcn.mat` (e.g. `QuadraticCost_Fcn.mat`) and another inheriting `Dynamics.mat`
(e.g. `ParafoilDynamics.mat`).


### File Descriptions

#### `SDDP.mat`

This file contains the function of the SDDP algorithm with back-tracking line search and regularization method 1 (please refer to the report).

#### `SDDP2.mat`

This file contains the function of the SDDP algorithm with back-tracking line search and regularization method 2 (please refer to the report).

#### `SDDP.mat`

This file contains the function of the SDDP algorithm with back-tracking line search and regularization method 3 (please refer to the report).

#### `backward_pass.mat`

This file contains a implementation of the backward pass used only in SDDP3.mat.

#### `SimulationModels/`

This directory contains the simulation models. Each model is saved as a separate file, which inherits from the abstract class Dynamics. Each model has its corresponding main file that runs the optimization and plots the results.

#### `HelperFunctions/`

This directory contains the helper functions used to generate the structure of the optimal solution and to simulate the dynamics with arbitrary control inputs.

## Usage

Before running the project, please ensure that all the necessary files are placed in the same folder. This includes the main script file, simulation files, and any additional helper functions.

## Acknowledgments

This project extends the code for classical DDP provided from the the following Git repository:

- [Repository A](https://github.com/maitreyakv/ddp-simulation.git): This repository served as a valuable reference for structuring the code.

Please note that this project is intended solely for educational purposes. I do not claim ownership of the original code, and any modifications or extensions made are done with the purpose of learning and experimentation.



## References

[1] E. Theodorou, Y. Tassa and E. Todorov, "Stochastic Differential Dynamic Programming," Proceedings of the 2010 American Control Conference, Baltimore, MD, USA, 2010, pp. 1125-1132, doi: 10.1109/ACC.2010.5530971.
