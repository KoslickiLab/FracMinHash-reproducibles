def num_expected_non_mutated_kmers(L, k, p):
    return L * (1-p)**k

for L in [10000, 100000, 1000000]:
    for k in [21, 51, 100]:
        for p in [0.001, 0.1, 0.2]:
            print(f'L: {L}, k:{k}, p:{p}, expetced num of nonmutated:{num_expected_non_mutated_kmers(L, k, p)}')

"""
$k = 21$  & 9792.1      & 1094.2    & 92.2    & 97920.9      & 10941.9     & 922.3    & 979208.7     & 109419.0    & 9223.3    \\ \hline
$k = 51$  & 9502.5      & 46.4      & 0.11     & 95025.4      & 463.8     & 1.1     & 950254.4     & 4638.4    & 11.4     \\ \hline
$k = 100$ & 9047.9      & 0.26     & 2.04E-6     & 90479.2      & 2.7      & 2.04E-5     & 904792.1     & 26.6    & 2.04E-4     \\ \hline
"""
