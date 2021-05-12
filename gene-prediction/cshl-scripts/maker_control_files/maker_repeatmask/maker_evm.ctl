#-----Transcript weights
evmtrans=10 #default weight for source unspecified est/alt_est alignments
evmtrans:blastn=0 #weight for blastn sourced alignments
evmtrans:est2genome=10 #weight for est2genome sourced alignments
evmtrans:tblastx=0 #weight for tblastx sourced alignments
evmtrans:cdna2genome=7 #weight for cdna2genome sourced alignments

#-----Protein weights
evmprot=10 #default weight for source unspecified protein alignments
evmprot:blastx=2 #weight for blastx sourced alignments
evmprot:protein2genome=10 #weight for protein2genome sourced alignments

#-----Abinitio Prediction weights
evmab=10 #default weight for source unspecified ab initio predictions
evmab:snap=10 #weight for snap sourced predictions
evmab:augustus=10 #weight for augustus sourced predictions
evmab:fgenesh=10 #weight for fgenesh sourced predictions
evmab:genemark=7 #weight for genemark sourced predictions
