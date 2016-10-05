 module load anaconda; snakemake --cluster "qsub {params.sge_opts}"  -p -s  ~/projects/HGSVG/hgsvg/stitching/Stitching.Snakefile  -j 22
