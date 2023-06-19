#!/usr/bin/python3.7
# coding=utf-8
"""
Validating the performance of PyF2F-Ruler by comparing
the obtained distance estimations published by Picco et al (Cell 2017)
with the distances estimations obtained running PyF2F-Ruler with the
same dataset.

We run a Chi-square test using the published set as the expected categoty (E)
and the set obtained with PyF2F-Ruler as the observed category (O).
"""
import sys
import numpy as np
import pandas as pd
import scipy.stats as stats

# Load Data for the C-ter
path_c = "../Cell_PyF2F_comparison/C_terminal.csv"
path_n = "../Cell_PyF2F_comparison/N_terminal.csv"

c_data_observed = pd.read_csv(path_c, usecols=["mu_py", "serr_mu_py"]).sort_values(['mu_py'])
n_data_observed = pd.read_csv(path_n, usecols=["mu_py", "serr_mu_py"]).sort_values(['mu_py'])
c_data_expected = pd.read_csv(path_c, usecols=["mu_cell", "serr_mu_cell"]).sort_values(['mu_cell'])
n_data_expected = pd.read_csv(path_n, usecols=["mu_cell", "serr_mu_cell"]).sort_values(['mu_cell'])

# Compare all the values (C & N ter)
# Observed values (PyF2F-Ruler distances) and Expected values (published dataset)
observed_values = np.concatenate((c_data_observed.mu_py.to_numpy(), n_data_observed.mu_py.to_numpy()), axis=0)
expected_values = np.concatenate((c_data_expected.mu_cell.to_numpy(), n_data_expected.mu_cell.to_numpy()), axis=0)

# Getting Frequencies from observed and expected values
# Binning data from min val (~13 nm) and max val (~36 nm)
obs_count, obs_bins = np.histogram(observed_values, bins=[13, 17, 21, 25, 29, 36])
exp_count, exp_bins = np.histogram(expected_values, bins=[13, 17, 21, 25, 29, 36])
chi_square_stat, p_value = stats.chisquare(obs_count, exp_count)

# Print Results
print(f"Chi-square test results in:\n"
      f"\n\tChi-square stat: {chi_square_stat}\n"
      f"\tp-value: {p_value}\n\n")

print("\nDone!\n")
sys.exit(0)

# END
