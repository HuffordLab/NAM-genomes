#!/bin/sh

SPLIT_PROG=/sonas-hs/ware/hpc/home/olson/split_gene_list.pl
FIXER_PROG=/sonas-hs/ware/hpc/home/olson/run_CDSfix_pipeline.pl

PROJ=$3
REGISTRY=$1
SPECIES=$2
ACCESSION=$3
CNT=$4

FIXROOT=/sonas-hs/ware/hpc_norepl/data/olson/NAM

if [[ ! $PROJ ]]
then
   echo "PROJ not defined, quit"
   exit
fi

if [[ ! $REGISTRY ]]
then
   echo "REGISTRY not defined, quit"
   exit
fi

if [[ ! $SPECIES ]]
then
   echo "SPECIES file not defined, quit"
   exit
fi

if [[ ! $ACCESSION ]]
then
   echo "ACCESSION file not defined, quit"
   exit
fi

if [[ ! $CNT ]]
then
   echo "Job Count not defined, quit"
   exit
fi

mkdir $PROJ
mkdir $PROJ/input
mkdir $PROJ/output
mkdir $PROJ/log



#split input file
$SPLIT_PROG -r $REGISTRY -s $SPECIES$ACCESSION -d $PROJ/input/ -c $CNT -l cshl_gene

echo "
#!/bin/sh

#$ -cwd
#$ -V
#$ -t 1-$CNT
#$ -N $PROJ 
#$ -o $PROJ/log/
#$ -e $PROJ/log/
#$ -l m_mem_free=4G 
#$ -l h_vmem=15G

echo "Task id is \$SGE_TASK_ID"

$FIXER_PROG --registry $REGISTRY --species $SPECIES$ACCESSION --guide $FIXROOT/fix_CDS/$ACCESSION.fixes_5.txt --stage 1 --genes $PROJ/input/$SPECIES$ACCESSION.\$SGE_TASK_ID
$FIXER_PROG --registry $REGISTRY --species $SPECIES$ACCESSION --guide $FIXROOT/fix_CDS/$ACCESSION.fixes_5.txt --stage 2 --genes $PROJ/input/$SPECIES$ACCESSION.\$SGE_TASK_ID
$FIXER_PROG --registry $REGISTRY --species $SPECIES$ACCESSION --guide $FIXROOT/fix_CDS/$ACCESSION.fixes_5.txt --stage 3 --genes $PROJ/input/$SPECIES$ACCESSION.\$SGE_TASK_ID
$FIXER_PROG --registry $REGISTRY --species $SPECIES$ACCESSION --guide $FIXROOT/fix_CDS/$ACCESSION.fixes_5.txt --stage 4 --genes $PROJ/input/$SPECIES$ACCESSION.\$SGE_TASK_ID
$FIXER_PROG --registry $REGISTRY --species $SPECIES$ACCESSION --guide $FIXROOT/fix_CDS/$ACCESSION.fixes_5.txt --stage 5 --genes $PROJ/input/$SPECIES$ACCESSION.\$SGE_TASK_ID

" >$PROJ.sh

chmod a+x $PROJ.sh

#submit job array
qsub $PROJ.sh
