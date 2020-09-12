import numpy as np
import pandas as pd
import random


"""
BUILD DATA FRAME OF PRIOR DRAWS FOR ABC
"""


#np.random.seed(214125)
np.random.seed(83242)

#n_draws = 3 * 20000
n_draws = 100000


#assume anc. Ne is 1000. median theta ~ 0.008. 4 * 1000 * 2e-6 = 0.008
#best estimate of per base mu is 3e-8. Going to use 2e-6
#best median estimate of recombination is r 1.6e-08 
#n = 3e-6/3e-8 = 67
#suggested rescaling is (1/2)*(1 - (1-2*r)^67) ~ 1.05e-6

#assume structural variant mutation rate is one of three possibilities
#mu_svs = [1e-10, 1e-9, 1e-8]
#mu_svs = [1e-9]

frac_0 = 0.1
s_binom = np.random.binomial(1, frac_0, n_draws).astype(bool)
sfs1_mean = -np.random.uniform(0, 0.1, n_draws)
sfs2_mean = -np.random.uniform(0, 0.1, n_draws)
sfs1_mean[s_binom] = 0
sfs2_mean[s_binom] = 0

Na = 1000
prior_df = pd.DataFrame({
    #fixed params
    'loci' : [20000000] * n_draws,
    'sample' : [26] * n_draws,
    'mu' : [2e-6] * n_draws,
    'rr' : [1.05e-6] * n_draws,
    'seeds' : np.random.randint(0, int(2**62) - 1, n_draws),
    'Na' : [Na] * n_draws, #ancestral pop size
    'N0' : np.exp(np.random.uniform(np.log(Na), np.log(20*Na), n_draws)).astype(int),
    'Nb' : np.random.randint(0.01*Na, Na, n_draws), #instant bottleneck pop size
    'B_t' : int(0.067*Na),  #time after bottleneck
    #'mu_sv' : np.random.choice(mu_svs, n_draws), #!!!!
    'mu_sv' : np.random.uniform(1e-10, 1e-7, n_draws), #!!!!
    'sfs1_mean' : sfs1_mean,
    'sfs1_shape' : np.random.uniform(0, 100, n_draws),
    'sfs2_mean' : sfs2_mean,
    'sfs2_shape' : np.random.uniform(0, 100, n_draws)
})

prior_df.to_csv('prior_df.csv')
