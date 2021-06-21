This is a small set of functions that help to compute different types of multivariate permutation entropies for time series analysis. As an example, the application to simulated multivariate fractional Brownian motion is provided.

1. [Getting started](#getting-started)
2. [Contributors and participation](#contribution-and-participation)
2. [Corresponding references](#corresponding-references)

Currently available:

+ Pooled Permutation Entropy ```pooled_permutation_entropy```
+ Multivariate Multiscale Permutation Entropy ```multivariate_multiscale_permutation_entropy```
+ Multivariate Weighted Permutation Entropy ```multivariate_weighted_permutation_entropy```
+ Multivariate Permutation Entropy based on PCA ```multivariate_permutation_entropy_pca```


# Getting started

## Create python virtual environment
1. install
```bash
pip install virtualenv
```

2. create an env (please change with local python version)
```bash
virtualenv ENV
```

3. Activate virtualenv by executing
```bash
. ENV/bin/activate
```

4. Close environment by executing 
```bash
deactivate
```

## Write/Freeze python dependency files
``` 
pip freeze > requirements.txt
``` 

## Install dependencies based on requirements file
```
pip install -r requirements.txt 
```

# Examples for Running 

The code below shows simple applications of ``multivariate_permutation_entropy.py``

.. code-block:: python

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


# Contributors and participation

* [Marisa Mohr](https://github.com/marisamohr)
* [Nils Finke](https://github.com/FinkeNils)


# Corresponding references
* https://github.com/nikdon/pyEntropy
* K.  Keller  and  H.  Lauffer,  “Symbolic  Analysis  of  High-DimensionalTime Series”, International Journal of Bifurcation and Chaos, vol. 13,no. 09, pp. 2657–2668, Sep. 2003.
* M.  Mohr,  F.  Wilhelm,  and  R.  Möller,  “On  the  Behaviour  of  Weighted Permutation  Entropy  on  Fractional  Brownian  Motion  in  the  Univariate and  Multivariate  Setting”, The  International  FLAIRS  Conference  Proceedings, vol. 34, Apr. 2021.
* F.  C.  Morabito,  D.  Labate,  F.  La  Foresta,  A.  Bramanti,  G.  Morabito,and  I.  Palamara,  “Multivariate  Multi-Scale  Permutation  Entropy  for Complexity  Analysis  of  Alzheimer’s  Disease  EEG ”, Entropy,  vol.  14,no. 7, pp. 1186–1202, Jul. 2012.
* M.  Mohr,  F.  Wilhelm,  M.  Hartwig,  R.  Möller,  and  K.  Keller,  “New Approaches  in  Ordinal  Pattern  Representations  for  Multivariate  TimeSeries,”  in Proceedings  of  the  33rd  International  Florida  Artificial Intelligence Research Society Conference (FLAIRS-33), 2020.
