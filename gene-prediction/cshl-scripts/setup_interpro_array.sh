#!/bin/sh

SPLIT_PROG=/sonas-hs/ware/hpc/home/weix/gramene-ensembl/scripts/misc-scripts/split-fasta.pl
INTERPROSCAN_PROG=/sonas-hs/ware/hpc/home/weix/my_interproscan/interproscan-5.38-76.0/interproscan.sh


PROJ=$1
INPUT=$2
CNT=$3
FMT=$4

if [[ ! $PROJ ]]
then
   echo "PROJ not defined, quit"
   exit
fi

if [[ ! $INPUT ]]
then
   echo "INPUT file not defined, quit"
   exit
fi

if [[ ! $CNT ]]
then
   echo "Job Count not defined, quit"
   exit
fi

if [[ ! $FMT ]]
then
  echo "using default FMT tsv"
  FMT="tsv"
fi

if [[  -d $PROJ ]]
then  
   echo "PROJ already exists, quit"
   exit
fi
 
if [[  ! -f $INPUT ]]
then  
   echo "INPUT not found, quit"
   exit
fi

   echo "split input $INPUT to $CNT jobs under $PROJ"

mkdir $PROJ
mkdir $PROJ/input
mkdir $PROJ/output
mkdir $PROJ/log


module load Java/11.0.2
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6


#split input file
perl $SPLIT_PROG -i $INPUT -d $PROJ/input/ -fn $CNT  -interproscan

cd $PROJ/

echo "
#!/bin/sh

#$ -cwd
#$ -V
#$ -t 1-$CNT
#$ -N $PROJ 
#$ -o test_log/
#$ -e test_log/
#$ -l m_mem_free=12G 
#$ -l h_vmem=15G

echo "Task id is \$SGE_TASK_ID"

$INTERPROSCAN_PROG -appl PfamA -iprlookup -f $FMT -i input/*.\$SGE_TASK_ID.fa -b output/$PROJ.\$SGE_TASK_ID
" >$PROJ.sh

chmod a+x $PROJ.sh

#submit job array
qsub $PROJ.sh
