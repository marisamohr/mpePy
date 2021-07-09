mpePy: A Python Package for Data Analysis with Multivariate Permutation Entropy for Time Series
===============================================================================================

``mpepy`` is a pure Python module  that implements data analysis methods based
on Bandt and Pompe's [#bandt_pompe]_ symbolic encoding scheme.

``mpepy`` implements the following data analysis methods:

- Pooled Permutation Entropy [#keller_lauffer]_; 
- Multivariate Multiscale Permutation Entropy [#morabito]_; 
- Multivariate Weighted Permutation Entropy [#mohr_a]_;
- Multivariate Ordinal Pattern Permutation Entropy [#mohr]_;
- Multivariate Permutation Entropy based on Principal Component Analysis [#mohr]_



Installing
==========

mpePy can be installed via the command line using

.. code-block:: console

   pip install mpepy

or you can directly clone its git repository:

.. code-block:: console

   git clone https://github.com/marisamohr/mpePy.git
   cd mpepy
   pip install -e .


Basic usage
===========


.. code-block:: python

    # Computing different multivariate permutation entropies for fractional Brownian motion.

    import subprocess
    import csv
    import pandas as pd
    import mpepy as mpe


    # Example of data simulation: multivariate fractional Brownian motion
    # usage of R-package 
    def simulateMultiFracBrownMotion(n, H_1, H_2, H_3, H_4, H_5, rho):
        output_file_name = './intermediate_output/MultiFracBrownMotionOutput.csv'
        subprocess.check_call(['Rscript', './intermediate_output/simulation_mfBm.R', str(n), str(H_1), str(H_2), str(H_3), str(H_4), str(H_5), str(rho), output_file_name], shell=False)
        arr = []
        with open(output_file_name, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                arr.append(row)
        mfbm = pd.DataFrame.from_records(arr)
        mfbm = mfbm.apply(pd.to_numeric)
        return mfbm
    # simulation
    mfbm = simulateMultiFracBrownMotion(2000, 0.3, 0.6, None, None, None, 0.0)
    mfbm = mfbm.T


    # Examples of multivariate permutation entropy calculation
    mpe.pooled_permutation_entropy(mfbm, order = 3 , delay = 1)
    mpe.multivariate_weighted_permutation_entropy(mfbm, order = 3 , delay = 1)
    mpe.multivariate_multiscale_permutation_entropy(mfbm, order = 3 , delay = 1, scale = 1)
    mpe.multivariate_ordinal_pattern_permutation_entropy(mfbm, order = 2 , delay = 1)
    mpe.multivariate_permutation_entropy_pca(mfbm, order = 2 , delay = 1, no_pc = 1)
    mpe.multivariate_permutation_entropy_pca(mfbm, order = 3 , delay = 5, no_pc = "all")


Contributors
============

- Marisa Mohr(https://github.com/marisamohr)
- Nils Finke(https://github.com/FinkeNils)



References
==========


.. [#bandt_pompe] Bandt, C., and Pompe, B. (2002). Permutation entropy: A Natural 
   Complexity Measure for Time Series. Physical Review Letters, 88, 174102.
.. [#mohr] Mohr, M., Wilhelm, F., Hartwig, M., Möller, R., and Keller, K. (2020). 
    New Approaches in Ordinal Pattern Representations for Multivariate Time Series. 
    In: Proceedings of the 33rd International Florida Artificial Intelligence 
    Research Society Conference (FLAIRS-33).
.. [#keller_lauffer] Keller, K., and Lauffer, H. (2003). Symbolic Analysis of 
    High-Dimensional Time Series. International Journal of Bifurcation and Chaos,
    vol. 13,no. 09, pp. 2657–2668.
.. [#mohr_a] Mohr, M., Wilhelm, F., and  Möller, R. (2021). On  the  Behaviour
    of Weighted Permutation Entropy on Fractional Brownian Motion in the Univariate
    and Multivariate Setting. The  International  FLAIRS  Conference  Proceedings,
    vol. 34.
.. [#morabito] Morabito, F.C., Labate, D., La  Foresta, F., Bramanti, A., Morabito, G.,
    and Palamara I. (2012). Multivariate  Multi-Scale  Permutation  Entropy  for 
    Complexity Analysis of Alzheimer’s Disease EEG. Entropy, vol. 14, no. 7.




