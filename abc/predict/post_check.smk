from fixed import *
import pandas as pd
import numpy as np


np.random.seed(125)

#can't run until completed by full pipeline because file changes so outfile names do too.
Smp = 20

abc_df = abc_df = pd.read_csv("data/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_df.csv", sep = ",").query('N0 > 26 and Nb > 2 and mu_sv > 0 and sfs1_mean > 0 and sfs2_mean > 0 and sfs1_shape > 0 and sfs2_shape > 0').groupby('window').apply(lambda x: x.sample(Smp) if len(x) > Smp else None).astype({'N0': 'int32', 'Nb': 'int32'})


#full_df = []
#mus = ["1e-10", "1e-09", "1e-08"]
#for mu_sv in mus: 
#    abc_df = pd.read_csv(f"data/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_fixedMu_{mu_sv}_df.csv", sep = ",").astype({'N0': 'int32', 'Nb': 'int32'}).query('N0 > 25 and Nb > 1 and sfs1_mean > 0 and sfs2_mean > 0 and sfs1_shape > 0 and sfs2_shape > 0').groupby('window').apply(lambda x: x.sample(Smp) if len(x) > Smp else None)   
#    print(f"mu_sv is {mu_sv}")
#    abc_df['mu_sv'] = mu_sv
#    full_df.append(abc_df)
#abc_df = pd.concat(full_df)


all_post = expand("data/postcheck_out/window__{window}_Na__1000_N0__{N0}_Nb__{Nb}_Bt__67_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt",
                        zip,
                        window = abc_df['window'],
                        N0 = abc_df['N0'],
                        Nb = abc_df['Nb'],
                        mu_sv = abc_df['mu_sv'],
                        sfs1_mean = abc_df['sfs1_mean'],
                        sfs1_shape = abc_df['sfs1_shape'],
                        sfs2_mean = abc_df['sfs2_mean'],
                        sfs2_shape = abc_df['sfs2_shape'],
                        )
#print("ALL POST", all_post)

sfs = "data/postcheck_out/post_sfs.txt"
param = "data/postcheck_out/post_params.txt"

rule all:
    input:
        all_post
        #sfs, param

rule post_check:
    input:
        posts = "data/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_df.csv",
        #posts = "data/NAM_abc_out/NAM_ABC_SFS-dup_inv_tra_ins_del_fixedMu_df.csv",
        slim ="nam_exons_rawsfs.slim"
    output:
       "data/postcheck_out/window__{window}_Na__{Na}_N0__{N0}_Nb__{Nb}_Bt__{Bt}_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt"
    params:
        window = "{window}",
        Na = "{Na}", N0 = "{N0}",Nb = "{Nb}", Bt = "{Bt}",
        mu_sv = "{mu_sv}",
        #mu_sv = mu_sv,
        sfs1shape="{sfs1_shape}", sfs1mean="{sfs1_mean}",
        sfs2shape="{sfs2_shape}", sfs2mean="{sfs2_mean}"
    shell:
        """
        slim  -define rr={rr} \
        -define n_size={sample} -define loci={loci} \
        -define mu={mu} -define mu_sv={params.mu_sv} \
        -define Na={params.Na} -define N0={params.N0} -define Nb={params.Nb} -define B_t={params.Bt} \
        -define sfs1_shape={params.sfs1shape} -define sfs1_mean={params.sfs1mean} \
        -define sfs2_shape={params.sfs2shape} -define sfs2_mean={params.sfs2mean} {input.slim} > {output}
        """


#rule sim_df:
#    input:
#        outs = expand("data/postcheck_out/window__{window}_Na__1000_N0__{N0}_Nb__{Nb}_Bt__67_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt",
#                        zip,
#                        window = abc_df['window'],
#                        N0 = abc_df['N0'],
#                        Nb = abc_df['Nb'],
#                        mu_sv = abc_df['mu_sv'],
#                        #mu_sv = mu_sv,
#                        sfs1_mean = abc_df['sfs1_mean'],
#                        sfs1_shape = abc_df['sfs1_shape'],
#                        sfs2_mean = abc_df['sfs2_mean'],
#                        sfs2_shape = abc_df['sfs2_shape'],
#                        )
#    output:
#        sfs = sfs,
#        param = param
#    run:
#        #shell("grep -h -B2 'SFS:' {input.outs} | grep 'SFS:' | cut -d ':' -f2- | sed 's/^ //g' > {output.sfs}"),
#        #shell("cat <( ls {input.outs} | sed 's/\_\_*/\t/g' | cut -f 2,4,6,8,10,12,14,16,18,20 | sort -u | sed 's/out\///g' ) <(ls {input.outs} | sed 's/\_\_*/\t/g' | cut -f 3,5,7,9,11,13,15,17,19,21) > {output.param}") 
#        
#        shell("> {output.sfs}")
#        shell("for i in `ls {input.outs}`; done grep -h -B2 'SFS:' $i | grep 'SFS:' | cut -d ':' -f2- | sed 's/^ //g' >> {output.sfs}; done"),
#        shell("cat <( ls {input.outs} | sed 's/\_\_*/\t/g' | cut -f 2,4,6,8,10,12,14,16,18,20 | sort -u | sed 's/out\///g' ) <(ls {input.outs} | sed 's/\_\_*/\t/g' | cut -f 3,5,7,9,11,13,15,17,19,21) > {output.param}") 
