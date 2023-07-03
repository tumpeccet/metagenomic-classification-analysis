import math
import numpy as np
from scipy.stats import wilcoxon
import csv
import argparse

parser = argparse.ArgumentParser(
    prog='Calculating statistics.',
    description='')
parser.add_argument("input", help="Input file path containing abundance results, must be in csv comma delimited format")
args = parser.parse_args()

def perform_permutation_test(observed_diff, permuted_diffs):
    observed_mean_diff = np.mean(observed_diff)
    permuted_mean_diffs = np.mean(permuted_diffs)
    print(observed_mean_diff)
    print(permuted_mean_diffs)
    p_value = np.mean(permuted_mean_diffs >= observed_mean_diff)
    return p_value

# Wilcoxon signed-rank test
def perform_wilcoxon_test(RCA_diff, RRCA_diff):
    test_statistic, p_value = wilcoxon(x=RCA_diff, y=RRCA_diff, alternative='greater')
    return test_statistic, p_value

with open(args.input, "r") as csv_abundances:
    # True results (TRA)
    TRA = []
    # kRRCA results
    kRRCA = []
    # kRCA results
    kRCA = []
    # mRRCA results
    mRRCA = []
    # mRCA results
    mRCA = []

    # skipping the headers
    csv_abundances.readline()

    # Split columns while reading
    for id, tax_id, name, tra, krrca, krca, mrca, mrrca in csv.reader(csv_abundances, delimiter=';'):
        TRA.append(float(tra))
        kRRCA.append(float(krrca))
        kRCA.append(float(krca))
        mRRCA.append(float(mrrca))
        mRCA.append(float(mrca))

    # True results (TRA)
    TRA = np.array(TRA)
    # kRRCA results
    kRRCA = np.array(kRRCA)
    # kRCA results
    kRCA = np.array(kRCA)
    # mRRCA results
    mRRCA = np.array(mRRCA)
    # mRCA results
    mRCA = np.array(mRCA)

    # Computing the differences for each algorithm
    kRCA_diff = TRA - kRCA
    kRRCA_diff = TRA - kRRCA

    mRCA_diff = TRA - mRCA
    mRRCA_diff = TRA - mRRCA

    print(kRCA_diff)
    print(kRRCA_diff)
    print(mRCA_diff)
    print(mRRCA_diff)

    # Perform permutation test
    num_permutations = int(math.pow(len(TRA), 2))
    print(num_permutations)
    permuted_diffs_kRRCA = np.zeros((num_permutations, len(TRA)))
    permuted_diffs_mRRCA = np.zeros((num_permutations, len(TRA)))

    for i in range(num_permutations):
        permuted_diffs_kRRCA[i] = np.random.permutation(kRRCA_diff - kRCA_diff)
        permuted_diffs_mRRCA[i] = np.random.permutation(mRRCA_diff - mRCA_diff)

    permutation_p_value_kRRCA = perform_permutation_test(kRRCA_diff - kRCA_diff, permuted_diffs_kRRCA)
    permutation_p_value_mRRCA = perform_permutation_test(mRRCA_diff - mRCA_diff, permuted_diffs_mRRCA)

    # Perform Wilcoxon signed-rank test
    wilcoxon_test_statistic_kRRCA, wilcoxon_p_value_kRRCA = perform_wilcoxon_test(kRRCA_diff, kRCA_diff)
    wilcoxon_test_statistic_mRRCA, wilcoxon_p_value_mRRCA = perform_wilcoxon_test(mRRCA_diff, mRCA_diff)

    # Print the test statistics and p-values
    print("kRCA vs kRRCA:")
    print("Test statistic:", wilcoxon_test_statistic_kRRCA)
    print("p-value (Wilcoxon):", wilcoxon_p_value_kRRCA)
    print("p-value (Permutation):", permutation_p_value_kRRCA)
    print()
    print("mRCA vs mRRCA:")
    print("Test statistic:", wilcoxon_test_statistic_mRRCA)
    print("p-value (Wilcoxon):", wilcoxon_p_value_mRRCA)
    print("p-value (Permutation):", permutation_p_value_mRRCA)


