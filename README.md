# Eco-evolution-in-3-guild-networks
The code and data presented in this repository constitutes the core for simulating the paper titled "Balancing selection pressures shapes trait evolution and resilience in three-guild networks."

## Code Description

- **Draw_Ternary.py**: Python code for drawing ternary diagrams.
- **coevolution.m**: Matlab code for calculates updated species traits based on coevolution.
- **LG_PHM.m**: Matlab function describing Community Dynamics (differential equations Eqn. 1a-1c).
- **get_FDQ.m**: Matlab code for computes functional diversity.
- **get_jacmat.m**: Matlab code for obtaining the Jacobian matrix.
- **get_matrix_rand.m**: Matlab code for generates a random interaction matrix. 
- **get_modularity.m**: Matlab code for calculates network modularity.
- **get_wNoda.m**: Matlab code for computes weighted nestedness (WNODA).
- **trait_match.m**: Matlab code for computes trait matching.
- **update_interaction.m**: Matlab code for updates interspecific interactions.

## Core Programs

- **Model_Core.m**: Matlab main program for generating data for ternary diagrams (requires saving corresponding network metric data as a txt file, e.g., “Directionality.txt”).
- **Draw_Ternary.py**: Python program used to read txt file data and draw it on ternary diagrams (note that "Draw_Ternary.py" and " Directionality.txt" must be placed in the same folder).

## Data Description

- **Measures of empirical ecological networks.csv**: Empirical ecological network data were obtained from the Web of Life database (www.web-of-life.es) and the Interaction Web Database (as referenced in [Andreazzi et al. 2017](https://doi.org/10.1086/692110)). The dataset includes 107 plant–pollinator mutualistic networks and 22 plant–herbivore antagonistic networks spanning a wide range of geographical regions.

## External Library

The Matlab library "BiMat" is utilized for analyzing network structures. It is publicly available [here](https://bimat.github.io/).
