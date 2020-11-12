#$ -S /bin/bash
#$ -cwd
#$ -l m_mem_free=4G
#$ -l tmp_free=300G
#$ -pe mpi 32

source /sonas-hs/ware/hpc/home/kchougul/.bash_profile
export LD_PRELOAD=/mnt/grid/ware/hpc/home/data/mcampbel/applications/openmpi-1.8.8/build/lib/libmpi.so
export AUGUSTUS_CONFIG_PATH=/mnt/grid/ware/hpc/home/data/mcampbel/applications/augustus-3.1/config

chromosome_fasta="$1" #unmasked genome fasta, setup repeatmasker options in maker_opts.ctl

/mnt/grid/ware/hpc/home/data/mcampbel/applications/openmpi-1.8.8/build/bin/mpiexec \
  -mca btl ^openib -n 32 \
  /mnt/grid/ware/hpc/home/data/mcampbel/applications/maker/bin/maker \
  -g $ $chromosome_fasta\
  -t 20 -fix_nucleotides -TMP $TMPDIR
