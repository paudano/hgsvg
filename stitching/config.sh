module unload python/3.5.2
module load python/2.7.3
source ~/environments/python27/bin/activate
export PATH=$HOME/software/bin:$PATH

module load pysam/0.8.4
module load numpy/1.11.0
module load biopython/1.63
module load pandas/0.20.3
module load gmp/6.1.2; module load mpfr/4.0.1; module load mpc/1.1.0;  module load gcc/latest
