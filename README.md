# Anticipating Critical Transitions using Early Warning Signals
This repository contains code to reproduce all analyses and figures in Dablander, Pichler, Cika, & Bacilieri ([2020](https://psyarxiv.com/5wc28)). Anticipating Critical Transitions in Psychological Systems using Early Warning Signals: Theoretical and Pratical Considerations.

The structure of this repository is as follows:

- *Code/* holds code relevant for the paper. Specifically:
    - *Code/simulate-data.R* has code to simulate from the model.
    - *Code/analyze-data.R* has code to analyze the simulation results.
    - *Code/figures.R* has code to create all figures. This is the main file.
    - *Code/helpers.R* collects a number of helper functions.
    - *Code/tests.R* has some tests for the code.
    - *Code/AUC.R* computes the Area under the Curve for the early warning signals.
    - *Code/boerlijst.R* has code for the replication and extension of Boerlijst et al. (2013).
- *Simulation-Results/* are available upon request due to the size of the files. Smaller simulations results are included.