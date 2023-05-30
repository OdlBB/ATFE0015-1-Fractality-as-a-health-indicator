# Fractality in Heart Rate Signals - Master's Thesis

This repository contains the code developed for my master's thesis on fractality in heart rate signals. The aim of the thesis is to explore the presence of fractal patterns in biological data and investigate their potential implications.

## Contents

- `ECG-data/`: This directory contains the first dataset used in the thesis. Due to privacy restrictions, the actual data files are not included, but information on how to obtain them is provided.

- `Heart-rates-diseases/`: This directory contains the second dataset used in the thesis. Due to privacy restrictions, the actual data files are not included, but information on how to obtain them is provided.

- `DFA/`: This directory contains notebooks performing various analysis on the datasets in order to identify the fractal pattern on a given signal. The main method used for this purpose is the DFA method.

- `Ivanov/`: This directory contains the notebooks used to build the Ivanov model, a stochastic model modeling the heart rate variability. Along with building the model, an analysis of the model and fitting the model to actual heart rate signals are performed in several notebooks.

- `GLM/`: This directory contains two notebooks that briefly overview the potential of using GLMs for modelling the heart rate.

- `FractalAnalysis.py`: this files contains the methods to identify the fractal pattern of a signal. Two methods are used:
    - The DFA method: it was implemented in three different ways:
        - `DFA`: The method was built from scratch using the various steps of the method
        - `DFA_fast`: The method is the exact same one but different tools were used to speed up the method
        - `DFA2`: The method was implemented in the `MDFA` package
    - The Power Spectrum Analysis is implemented in the `PowerSpectrumAnalysis` method

- `FractalModel.py`: provides the method `IvanovModel` which returns a simulated signal of the heart rate.


## References

- P.C. Ivanov et al. “Stochastic feedback and the regulation of biological rhythms”. In: Europhysics letters 43 (1998), pp. 363–368.
- C.-K. Peng et al. “Mosaic organization of DNA nucleotides”. In: PHYSICAL REVIEW E 42.2 (1994), pp. 1685–1689.
- Marlene Müller. “Generalized linear model”. In: Handbook of Computational Statistics (2004), pp. 681–709.
