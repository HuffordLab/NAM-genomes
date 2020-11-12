import numpy as np
import pandas as pd
import random

#TO TRIGGER NEW RUN OF PIPELINE, RUN:
#$ python src/prior_df.py
#NOTE THERE IS A RANDOM SEED SET IN THE PRIOR_DF.PY SCRIPT

#assume anc. Ne is 1000. median theta ~ 0.008. 4 * 1000 * 2e-6 = 0.008
#best estimate of per base mu is 3e-8. Going to use 2e-6
#best median estimate of recombination is r 1.6e-08 
#n = 3e-6/3e-8 = 67
#suggested rescaling is (1/2)*(1 - (1-2*r)^67) ~ 1.05e-6

#fixed params
#loci = 20000000 #!!!
#sample = 26
#mu = 2e-6
#rr = 1.05e-6

#priors
#n_draws = 150000
#n_draws = 90000

#random.seed(214125) #random seed per sim is reported, so no reason to make this stage reproducible afaikt
#seeds = np.random.randint(0, int(2**62) - 1, n_draws)
#low informative prior set
#Na = 1000 #ancestral pop size
#N0 = np.random.randint(Na, 10*Na, n_draws) #modern pop size
#N0 = np.exp(np.random.uniform(np.log(Na), np.log(20*Na), n_draws)).astype(int)
#Nb = np.random.randint(0.05*Na, Na, n_draws) #instant bottleneck pop size
#B_t = int(0.067*Na)  #time after bottleneck
#mu_sv = np.random.uniform(0, 5e-8, n_draws)#.astype(np.half)
#sfs1_mean = -np.random.uniform(0, 0.1, n_draws)#.astype(np.half) #DFE negative only!
#sfs1_shape = np.random.uniform(2, 100, n_draws)#.astype(np.half) #DFE positive only!
#sfs2_mean = -np.random.uniform(0, 0.1, n_draws)#.astype(np.half) #DFE negative only!
#sfs2_shape = np.random.uniform(2, 100, n_draws)#.astype(np.half) #DFE positive only!

prior_df = pd.read_csv('prior_df.csv').astype({'Na': 'int64', 'N0': 'int64', 'Nb': 'int64'})

stats_out = [ f"abc_out/seed__{int(row['seeds'])}_Na__{int(row['Na'])}_N0__{int(row['N0'])}_Nb__{int(row['Nb'])}_Bt__{int(row['B_t'])}_musv__{row['mu_sv']}_sfs1shape__{row['sfs1_shape']}_sfs1mean__{row['sfs1_mean']}_sfs2shape__{row['sfs2_shape']}_sfs2mean__{row['sfs2_mean']}_sumstats.txt" for index, row in prior_df.iterrows()]


#fixed params
loci = list(prior_df['loci'])[0]
sample =  list(prior_df['sample'])[0]
mu =  list(prior_df['mu'])[0]
rr =  list(prior_df['rr'])[0]


#for i in range(n_draws):
#    param_str = f"abc_out/seed__{seeds[i]}_Na__{Na}_N0__{N0[i]}_Nb__{Nb[i]}_Bt__{B_t}_musv__{mu_sv[i]}_sfs1shape__{sfs1_shape[i]}_sfs1mean__{sfs1_mean[i]}_sfs2shape__{sfs2_shape[i]}_sfs2mean__{sfs2_mean[i]}_sumstats.txt"
#    stats_out.append(param_str)

rule all:
    input:
        stats_out
        #"pop_sfs.txt",
        #"pop.txt"

rule run_slim:
    input:
        slim = "src/nam_exons_rawsfs.slim",
        df = "src/prior_df.py"
    output:
        slim_out = "abc_out/seed__{seeds}_Na__{Na}_N0__{N0}_Nb__{Nb}_Bt__{B_t}_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt",
    params:
        seeds = "{seeds}" ,
        Na = "{Na}", N0 = "{N0}",Nb = "{Nb}", Bt = "{B_t}",
        mu_sv = "{mu_sv}",
        sfs1shape="{sfs1_shape}", sfs1mean="{sfs1_mean}", 
        sfs2shape="{sfs2_shape}", sfs2mean="{sfs2_mean}" 
        #pneutral="{p_neutral}", psfs1="{p_sfs1}"
    shell:
        """
        slim -s {params.seeds} -define rr={rr} \
        -define n_size={sample} -define loci={loci} \
        -define mu={mu} -define mu_sv={params.mu_sv} \
        -define Na={params.Na} -define N0={params.N0} -define Nb={params.Nb} -define B_t={params.Bt} \
        -define sfs1_shape={params.sfs1shape} -define sfs1_mean={params.sfs1mean} \
        -define sfs2_shape={params.sfs2shape} -define sfs2_mean={params.sfs2mean} {input.slim} > {output.slim_out}
        """


#for i in abc_out/*Na__1000_*; do grep -h -B2 "SFS:" $i | grep "SFS:" | cut -d ":" -f2- | sed 's/^ //g'; done > pop_sfs.txt &
#for i in abc_out/*Na__1000_*; do grep -h -B2 "SFS:" $i | grep -h -A1 "#Na" | grep -v "\-\-" | grep -v "#Na"; done > pop.txt

rule sim_df:
    input:
        expand("abc_out/seed__{seeds}_Na__{Na}_N0__{N0}_Nb__{Nb}_Bt__{B_t}_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt", 
                        zip,
                        seeds =  prior_df['seeds'],
                        Na = prior_df['Na'],
                        N0 = prior_df['N0'],
                        Nb = prior_df['Nb'],
                        B_t = prior_df['B_t'],
                        mu_sv = prior_df['mu_sv'],
                        sfs1_mean = prior_df['sfs1_mean'],
                        sfs1_shape = prior_df['sfs1_shape'],
                        sfs2_mean = prior_df['sfs2_mean'],
                        sfs2_shape = prior_df['sfs2_shape'],
                        )
    output:
        sfs = "pop_sfs.txt",
        param = "pop.txt"
    shell:
        """
            grep -h -B2 "SFS:" {input} | grep "SFS:" | cut -d ":" -f2- | sed 's/^ //g' > {output.sfs}
            grep -h -B2 "SFS:" {input} | grep -h -A1 "#Na" | grep -v "\-\-" | grep -v "#Na" > {output.param} 
        """

