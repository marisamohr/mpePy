import itertools
import numpy as np
import pandas as pd
import permutation_entropy as ent
import copy
from sklearn.decomposition import PCA


""" 
Calculation of different multivariate permutation entropy, e.g.,
- pooled permutation entropy,
- multivariate weighted permutation entropy, 
- multivariate multiscale permutation entropy,
- multivariate ordinal pattern permutation entropy, and
- multivariate permutation entropy based on principal component analysis 

""" 

def pooled_permutation_entropy(mts, order, delay):
    """ 
    based on K.  Keller  and  H.  Lauffer,  “Symbolic  Analysis  of  High-DimensionalTime Series”, International Journal of Bifurcation and Chaos, vol. 13,no. 09, pp. 2657–2668, Sep. 2003.
    """ 
    no_of_dim_m = mts.shape[0]
    no_of_samples_T = mts.shape[1]
    permutations = np.array(list(itertools.permutations(range(order))))

    total_no_of_comparisons_delta = (no_of_samples_T - (order * delay - delay)) * no_of_dim_m
    start_pattern = order * delay - delay + 1
    freq_matrix = []
    total_freq_matrix = []
    total_pattern_arr = []
    total_pattern_arr_including_values = []

    # iterate over each row in matrix/each dimension
    for idx_dim in range(0, no_of_dim_m):
        series = mts.iloc[idx_dim][:]
        series = series.values
        
        # create pattern array for time series
        pattern_arr = []
        for idx_sample in range(start_pattern, no_of_samples_T):
            extract = series[idx_sample-(start_pattern):idx_sample:delay]
            extract_sorted = copy.copy(extract)
            extract_sorted.sort() 

            pattern = []

            # do index lookup to reduce from value representation (e.g. 1 7 5) to "pattern" represenation (e.g. 0 2 1)
            for number in extract_sorted:
                idx_in_extract = np.where(extract == number)
                idx_in_extract = idx_in_extract[0][0]

                pattern.append(idx_in_extract)

            pattern_arr.append(pattern)
            total_pattern_arr.append(pattern)

            # calc varinance for pattern and store
            total_pattern_arr_including_values.append({
                'pattern': pattern,
                'values': extract
            })
        
        # case 1: count different patterns per permutation per dimension
        for permuation in permutations:
            count_local = 0

            for pattern in pattern_arr:
                if np.array_equal(pattern,permuation):
                    count_local = count_local + 1

            # determine counts for each permutation
            pattern_dict = {
                'timeseries_idx': idx_dim,
                'permuation': permuation,
                'count': count_local,
                'relative_freq': count_local / total_no_of_comparisons_delta
            }
            freq_matrix.append(pattern_dict)

    # case 2: count different patterns per permutation independet of the dimension
    for permuation in permutations:
        count_global = 0

        for pattern_dict in total_pattern_arr_including_values:
            pattern = pattern_dict['pattern']
            if np.array_equal(pattern,permuation):
                count_global = count_global + 1

        total_pattern_dict = {
            'permuation': permuation,
            'count': count_global,
            'marginal_relative_freq': count_global / total_no_of_comparisons_delta
        } 
        total_freq_matrix.append(total_pattern_dict)

    ppe = 0 
    for marginal_dict in total_freq_matrix:
        freq = marginal_dict['marginal_relative_freq']
        ppe = ppe + (freq * np.log(freq))
    ppe = ppe * -1

    return ppe



def multivariate_weighted_permutation_entropy(mts, order, delay):
    """ 
    based on M.  Mohr,  F.  Wilhelm,  and  R.  Möller,  “On  the  Behaviour  of  Weighted Permutation  Entropy  on  Fractional  Brownian  Motion  in  the  Univariate and  Multivariate  Setting”, The  International  FLAIRS  Conference  Proceedings, vol. 34, Apr. 2021.
    """ 
    no_of_dim_m = mts.shape[0]
    no_of_samples_T = mts.shape[1]
    permutations = np.array(list(itertools.permutations(range(order))))

    total_no_of_comparisons_delta = (no_of_samples_T - (order * delay - delay)) * no_of_dim_m
    start_pattern = order * delay - delay + 1
    freq_matrix = []
    total_freq_matrix = []
    total_pattern_arr = []
    total_pattern_arr_including_values = []

    # iterate over each row in matrix/each dimension
    for idx_dim in range(0, no_of_dim_m):
        series = mts.iloc[idx_dim][:]
        series = series.values
        
        # create pattern array for time series
        pattern_arr = []
        for idx_sample in range(start_pattern, no_of_samples_T):
            extract = series[idx_sample-(start_pattern):idx_sample:delay]
            extract_sorted = copy.copy(extract)
            extract_sorted.sort() 

            pattern = []

            # do index lookup to reduce from value representation (e.g. 1 7 5) to "pattern" represenation (e.g. 0 2 1)
            for number in extract_sorted:
                idx_in_extract = np.where(extract == number)
                idx_in_extract = idx_in_extract[0][0]

                pattern.append(idx_in_extract)

            pattern_arr.append(pattern)
            total_pattern_arr.append(pattern)

            # calc varinance for pattern and store
            variance = np.var(extract)
            total_pattern_arr_including_values.append({
                'pattern': pattern,
                'values': extract,
                'variance': variance
            })

        
        # case 1: count different patterns per permutation per dimension
        for permuation in permutations:
            count_local = 0

            for pattern in pattern_arr:
                if np.array_equal(pattern,permuation):
                    count_local = count_local + 1

            # determine counts for each permutation
            pattern_dict = {
                'timeseries_idx': idx_dim,
                'permuation': permuation,
                'count': count_local,
                'relative_freq': count_local / total_no_of_comparisons_delta
            }
            freq_matrix.append(pattern_dict)

    # case 2: count different patterns per permutation independet of the dimension
    for permuation in permutations:
        count_global = 0
        sum_variance = 0

        for pattern_dict in total_pattern_arr_including_values:
            pattern = pattern_dict['pattern']
            if np.array_equal(pattern,permuation):
                sum_variance += pattern_dict['variance']
                count_global = count_global + 1

        total_pattern_dict = {
            'permuation': permuation,
            'count': count_global,
            'sum_variance': sum_variance,
            'marginal_relative_freq': count_global / total_no_of_comparisons_delta
        } 
        total_freq_matrix.append(total_pattern_dict)

    total_variance = 0
    for marginal_dict in total_freq_matrix:
        sum_variance_per_perm = marginal_dict['sum_variance']
        total_variance += sum_variance_per_perm

    mwpe = 0
    for marginal_dict in total_freq_matrix:
        sum_variance_per_perm = marginal_dict['sum_variance']
        freq_var = sum_variance_per_perm / total_variance
        mwpe = mwpe + (freq_var * np.log(freq_var))


    mwpe = mwpe * -1

    return  mwpe




def multivariate_multiscale_permutation_entropy(mts, order, delay, scale):
    """ 
    based on F.  C.  Morabito,  D.  Labate,  F.  La  Foresta,  A.  Bramanti,  G.  Morabito,and  I.  Palamara,  “Multivariate  Multi-Scale  Permutation  Entropy  for Complexity  Analysis  of  Alzheimer’s  Disease  EEG ”, Entropy,  vol.  14,no. 7, pp. 1186–1202, Jul. 2012.
    """ 
    res_scaled = pd.DataFrame()
    ts_length = mts.shape[0]
    iter_length = int(ts_length/scale)
    for i in range(0, iter_length):
        row_at_i = mts.iloc[i*scale:((i*scale)+scale)]
        mean_row_at_i = row_at_i.mean()
        res_scaled = res_scaled.append(mean_row_at_i, ignore_index = True)

    mmspe = pooled_permutation_entropy(res_scaled, order , delay)
    
    return  mmspe




def multivariate_ordinal_pattern_permutation_entropy(mts, order, delay):
    """ 
    based on M.  Mohr,  F.  Wilhelm,  M.  Hartwig,  R.  Möller,  and  K.  Keller,  “New Approaches  in  Ordinal  Pattern  Representations  for  Multivariate  TimeSeries,”  in Proceedings  of  the  33rd  International  Florida  Artificial Intelligence Research Society Conference (FLAIRS-33), 2020
    """ 
    window_size = (order - 1) * delay + 1
    no_of_dim_m = mts.shape[0]
    no_of_samples_T = mts.shape[1]
    #permutations = np.array(list(itertools.permutations(range(order))))

    multivariate_pattern_counts_as_dict = {}

    # iterative window-wise through multivariate time series
    for t in range(window_size, no_of_samples_T, 1):
        # determine column index to extract from time series
        window_column_arr = np.arange((t-window_size), t)
        window_column_arr_skipped_delay = []
        delay_counter = 0
        for column_idx in window_column_arr:
            if delay_counter % delay == 0:
                window_column_arr_skipped_delay.append(column_idx)
            delay_counter += 1

        window = mts[window_column_arr_skipped_delay]

        # determine ordinal pattern per dimension
        multivariate_pattern = []
        for dim in range(0, no_of_dim_m):

            # determine pattern for window
            pattern = []
            dim_window = window.loc[dim, :].to_numpy()
            dim_window_sorted = copy.copy(dim_window)
            dim_window_sorted.sort() 

            # do index lookup to reduce from value representation (e.g. 1 7 5) to "pattern" represenation (e.g. 0 2 1)
            for number in dim_window_sorted:
                idx_in_window = np.where(dim_window == number)
                idx_in_window = idx_in_window[0][0]
                pattern.append(idx_in_window)

            multivariate_pattern.append(pattern)

        # determine unique indicator for multivariate pattern
        hashmult = np.power(order, np.arange((order*no_of_dim_m)))
        hashmult_matrix = hashmult.reshape((no_of_dim_m, order))
        hashval = (np.multiply(multivariate_pattern, hashmult_matrix)).sum()

        # append hashval as indicatior for multivariate pattern to counter
        contains_pattern = str(hashval) in multivariate_pattern_counts_as_dict.keys()
        if contains_pattern:
            multivariate_pattern_counts_as_dict[str(hashval)] = multivariate_pattern_counts_as_dict[str(hashval)] + 1
        else:
            multivariate_pattern_counts_as_dict[str(hashval)] = 1

    # determine frequency array based on patterns
    freq_arr = []
    for key in multivariate_pattern_counts_as_dict.keys():
        m_pattern_count = multivariate_pattern_counts_as_dict[key]
        freq = (m_pattern_count / (no_of_samples_T - delay*(order-1)))
        freq_arr.append(freq)

    # calculate moppe based on frequency array
    moppe = 0
    for freq in freq_arr:
        moppe += freq * np.log(freq)
    moppe = moppe * -1
    
    return moppe




def multivariate_permutation_entropy_pca(mts, order, delay, no_pc):
    """ 
    based on M.  Mohr,  F.  Wilhelm,  M.  Hartwig,  R.  Möller,  and  K.  Keller,  “New Approaches  in  Ordinal  Pattern  Representations  for  Multivariate  TimeSeries,”  in Proceedings  of  the  33rd  International  Florida  Artificial Intelligence Research Society Conference (FLAIRS-33), 2020
    """ 
    dimensions = mts.shape[0]
    # perform principal component analysis
    pca = PCA(n_components=dimensions)
    mts_pca_arr = pca.fit_transform(mts.T)    
    mts_pca = pd.DataFrame(mts_pca_arr)  
    # calculate permutation entropy based on principal component analysis
    if no_pc == "all":
        # calculate mpe-ppe, i.e., pooled permutation entropy based on all principal components
        mpe_pca = pooled_permutation_entropy(mts_pca.T, order , delay)
    else:
        # calculate mpe-pca_i, i.e., permutation entropy based on the i-th principal component
        pc = mts_pca.loc[:,no_pc].to_numpy()
        mpe_pca = ent.permutation_entropy(pc, order , delay )
    
    return  mpe_pca

