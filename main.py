import subprocess
import csv
import pandas as pd
import multivariate_permutation_entropy as mpe


# Example of data simulation: multivariate fractional Brownian motion
# usage of R-package 
def simulateMultiFracBrownMotion(n, H_1, H_2, H_3, H_4, H_5, rho):
    output_file_name = './intermediate_output/MultiFracBrownMotionOutput.csv'
    subprocess.check_call(['Rscript', './simulation_mfBm.R', str(n), str(H_1), str(H_2), str(H_3), str(H_4), str(H_5), str(rho), output_file_name], shell=False)
    arr = []
    with open(output_file_name, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            arr.append(row)
    mfbm = pd.DataFrame.from_records(arr)
    mfbm = mfbm.apply(pd.to_numeric)
    return mfbm
# simulation
mfbm = simulateMultiFracBrownMotion(1000, 0.3, 0.6, None, None, None, 0.0)


# Examples of multivariate permutation entropy calculation
mpe.pooled_permutation_entropy(mfbm.T, order = 3 , delay = 1)
mpe.multivariate_weighted_permutation_entropy(mfbm.T, order = 3 , delay = 1)
mpe.multivariate_multiscale_permutation_entropy(mfbm.T, order = 3 , delay = 1, scale = 1)
mpe.multivariate_permutation_entropy_pca(mfbm.T, order = 3 , delay = 1, no_pc = 1)
mpe.multivariate_permutation_entropy_pca(mfbm.T, order = 3 , delay = 5, no_pc = "all")


