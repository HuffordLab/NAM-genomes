#$ -S /bin/bash
#$ -cwd
#$ -pe threads 16 -l m_mem_free=3.8G
#$ -l tmp_free=300G

source /sonas-hs/it/hpc/home/easybuild/lmod-setup.sh
module load Anaconda3/5.3.0
source activate mypasa

mikadao_combined_fasta="$1"
genome_fasta="$2" #repeatmasked
alignAssembly_config="$3"
acc_list="$4"

filename=$(basename $mikadao_combined_fasta)
cln=${filename}.cln
clean=${filename}.clean

$PASAHOME/bin/seqclean $mikadao_combined_fasta \
  -r $cln \
  -o $clean
  
$PASAHOME/Launch_PASA_pipeline.pl \
  -C -R -T \
  -c $alignAssembly_config \
  -g $genome_fasta \ 
  -t $clean  
  -u $mikado_combined_fasta \
  -f $acc_list \
  --ALIGNERS blat,gmap \
  --CPU 16
