#$ -S /bin/bash
#$ -cwd
#$ -N luj-maize.flc.iso.est.B73_AB10_mikado.combined.fasta
#$ -o /mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B73_AB10/pasa_run/output/maize.flc.iso.est.B73_AB10_mikado.combined.fasta-out.txt
#$ -e /mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B73_AB10/pasa_run/error/maize.flc.iso.est.B73_AB10_mikado.combined.fasta-err.txt
#$ -terse
###$ -l virtual_free=3.8G
#$ -l m_mem_free=3.8G
###$ –l h_vmem=3.8G



echo '=================================================='
# print local SGE vars
echo JOB_ID=$JOB_ID
echo QUEUE=$QUEUE
echo SGE_TASK_ID=$SGE_TASK_ID
echo TMPDIR=$TMPDIR
echo PWD=$PWD
echo SUBMIT_TIME=`date`
echo '=================================================='

#$ -v PATH,CLASSPATH,JARPATH,PERL5LIB,LD_LIBRARY_PATH,DATANR_HOME,DATA_HOME,DATAFC_HOME



# staging


# start of template

# specify that the job requires 3GB of memory for each process
#$ -pe threads 16 -l m_mem_free=3.8G
# make sure that there is at least 20G of tmp space on the node
#$ -l tmp_free=300G

source /sonas-hs/it/hpc/home/easybuild/lmod-setup.sh
module load Anaconda3/5.3.0
source activate mypasa 

$PASAHOME/bin/seqclean /mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B73_AB10/maize.flc.iso.est.B73_AB10_mikado.combined.fasta -r maize.flc.iso.est.B73_AB10_mikado.combined.fasta.cln -o maize.flc.iso.est.B73_AB10_mikado.combined.fasta.clean || { echo JOB_STATUS=ERROR 1>&2; exit; }
$PASAHOME/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g /mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/NAM_Canu1.8_new_runs/B73Ab10/maker_final_annotations/B73Ab10.maker.repeatmasked.fasta -t maize.flc.iso.est.B73_AB10_mikado.combined.fasta.clean -T -u /mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B73_AB10/maize.flc.iso.est.B73_AB10_mikado.combined.fasta -f /mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/FL.acc.list --ALIGNERS blat,gmap --CPU 16




# end of template
echo '=================================================='
echo END_TIME=`date`
echo '=================================================='

#################################################



















