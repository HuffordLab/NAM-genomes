#$ -S /bin/bash
#$ -cwd
#$ -pe threads 16 -l m_mem_free=3.8G
#$ -l tmp_free=300G

source /sonas-hs/it/hpc/home/easybuild/lmod-setup.sh
module load Anaconda3/5.3.0
source activate mypasa

pasa_rnd1_gff3="$1" #gff3 from first round pasa run
genome_fasta="$2" #repeatmasked
clean="$3" #cleaned sequence from step 1
alignAssembly_config="$4"
annotCompare_config="$5"

$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi \
  -c $alignAssembly_config \
  -g $genome_fasta \
  -P $pasa_rnd1_gff3

$PASAHOME/Launch_PASA_pipeline.pl \
  -A --CPU 16 \
  -c $annotCompare_config \
  -g $genome_fasta \
  -t $clean
