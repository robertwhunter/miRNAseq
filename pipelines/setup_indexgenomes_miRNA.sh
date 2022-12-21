#!/bin/sh
# index genomes (only once)

# Set grid Engine options:
#$ -N index_genomes
#$ -wd <insert wd>
#$ -e <insert wd>/logs/
#$ -o <insert wd>/logs/
#$ -l h_rt=01:00:00
#$ -l h_vmem=12G
#$ -pe sharedmem 2
#$ -m bea -M <insert email>
#$ -V

set_conda_dir=$DIR_miniconda
source $set_conda_dir/bin/activate env_MAPmiR

# U>T in miRBase files
sed 's/U/T/g' miRNA_U.fa > miRNA.fa
sed 's/U/T/g' stemloop_U.fa > stemloop.fa
rm *_U.fa

# build indeces
parallel "bowtie-build {} {.}" ::: *.fa

conda deactivate
