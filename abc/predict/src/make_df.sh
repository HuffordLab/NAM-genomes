#data/sfs/NAM_SFS-DUP_SV.txt
#data/sfs/NAM_SFS_INS_DEL-gte2Kb_SV.txt
#data/sfs/NAM_SFS_INS_DEL-lt2Kb_SV.txt
#data/sfs/NAM_SFS-INV_SV.txt
#data/sfs/NAM_SFS-TRA_SV.txt

#paste data/sfs/NAM_SFS_4fold_vcf.txt  data/sfs/NAM_SFS_0fold_vcf.txt  data/sfs/NAM_SFS_SV.txt > data/sfs/NAM_all.txt
paste data/sfs/NAM_SFS_4fold_vcf.txt  data/sfs/NAM_SFS_0fold_vcf.txt  data/sfs/NAM_SFS-DUP_SV.txt > data/sfs/NAM_all_DUP.txt

paste data/sfs/NAM_SFS_4fold_vcf.txt  data/sfs/NAM_SFS_0fold_vcf.txt  data/sfs/NAM_SFS_INS_DEL-gte2Kb_SV.txt > data/sfs/NAM_all_INDEL_BIG.txt

paste data/sfs/NAM_SFS_4fold_vcf.txt  data/sfs/NAM_SFS_0fold_vcf.txt  data/sfs/NAM_SFS_INS_DEL-lt2Kb_SV.txt > data/sfs/NAM_all_INDEL_SMALL.txt

paste data/sfs/NAM_SFS_4fold_vcf.txt  data/sfs/NAM_SFS_0fold_vcf.txt  data/sfs/NAM_SFS-INV_SV.txt > data/sfs/NAM_all_INV.txt

paste data/sfs/NAM_SFS_4fold_vcf.txt  data/sfs/NAM_SFS_0fold_vcf.txt  data/sfs/NAM_SFS-TRA_SV.txt > data/sfs/NAM_all_TRA.txt


