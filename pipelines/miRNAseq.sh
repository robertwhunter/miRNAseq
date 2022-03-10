#!/bin/sh

#$ -N <insert job name>_miRNAseq
#$ -wd <insert wd>
#$ -e <insert wd>/logs/
#$ -o <insert wd>/logs/
#$ -l h_rt=08:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 6
#$ -m bea -M <insert email>
#$ -V


####################################################################
#### READ-ME #######################################################
####################################################################
##
## Input files should be trimmed .fq files in `trim` directory
## Genome directory must contain .fa and Bowtie index files
##
## Before running, set grid engine options and script parameters
## NB Conda directory and genome directory are set here as environmental variables
##
## Ensure set_cores matches no. of cores in Grid Engine
##


####################################################################
#### SET-UP ########################################################
####################################################################
##
##

set_cores=6
set_dummyrun=1
set_dummysize=400000
set_conda_dir=$DIR_miniconda    # no trailing slash
set_genomes_dir=$DIR_genome_miR # no trailing slash
mod1=1 # QC
mod2=1 # trim
mod3=1 # map 

#### Initialise the environment modules

. /etc/profile.d/modules.sh


#### Write out

echo "Starting miRNAseq mapping pipeline"
echo "" 
echo ""
echo "Running modules:" 
echo "----------------"
if (( $mod1 == 1 )); then echo "Module 1 - QC"; fi
if (( $mod2 == 1 )); then echo "Module 2 - trim"; fi
if (( $mod3 == 1 )); then echo "Module 3 - map"; fi
echo "" 
echo ""
echo "Script parameters:"
echo "----------------"
echo "Genome dir = $set_genomes_dir"
echo ""
echo ""
echo "Input fastq files:"
echo "----------------"
ls -ho *.fastq
echo ""
echo ""
if (( $set_dummyrun == 1 )); then echo "Running as dummy run"; fi


####################################################################
#### COPY FASTQ FILES ##############################################
####################################################################
##
##

t="$(date)" 
echo "...starting fastq prep at $t..."

mkdir fastq/

if (( $set_dummyrun == 1 )) 

then
  parallel "head -$set_dummysize {} > fastq/{}" ::: *.fastq 

else
  parallel "cp {} fastq/{}" ::: *.fastq 

fi

source $set_conda_dir/bin/activate env_UMI

#### Module 1
if (( $mod1 == 1 ))
then

#### QC on raw files

mkdir raw_qc
mv fastq/*.fastq raw_qc

parallel "fastqc {}" ::: raw_qc/*.fastq
mv raw_qc/*.fastq fastq

multiqc raw_qc -n raw_multiqc
mkdir raw_multiqc
mv raw_multiqc.html raw_multiqc_data raw_multiqc

fi

####################################################################
#### TRIM ##########################################################
####################################################################
##
##

#### Module 2
if (( $mod2 == 1 ))
then

t="$(date)" 
echo "...starting trimming at $t..."

source $set_conda_dir/bin/activate env_UMI

#### UMI processing

for f in fastq/*.fastq; do

umi_tools extract \
  --extract-method=regex \
  --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)' \
  --stdin=${f} \
  --log=${f%.*}.umi.log \
  --stdout=${f%.*}.umi 

done

#### Trim with cutadapt

parallel \
  "cutadapt \
  -a AACTGTAGGCACCATCAAT \
  --minimum-length 18 --maximum-length 30 \
  -o {.}.trim.fq \
  {} \
  > {.}.cutadapt.log" ::: fastq/*.umi

mv fastq/*log logs

  if (( $mod1 == 1 ))
  then

  #### QC on trimmed files

  mkdir trim_qc
  mv fastq/*.trim.fq trim_qc

  parallel "fastqc {}" ::: trim_qc/*.trim.fq
  mkdir trim
  mv trim_qc/*.trim.fq trim

  multiqc trim_qc -n trim_multiqc
  mkdir trim_multiqc
  mv trim_multiqc.html trim_multiqc_data trim_multiqc
  
  else
  
  mkdir trim
  mv fastq/*.trim.fq trim
  
  fi
  
#### Get final library sizes

touch trim/lib_size.txt
echo "library reads" > trim/lib_size.txt

for f in trim/*.fq; do
  numlines=$(wc -l < $f)
  numreads=$(expr $numlines / 4)
  echo "${f##*/} ${numreads}" >> trim/lib_size.txt
done  

fi

conda deactivate



####################################################################
#### MAP ###########################################################
####################################################################
##
##

#### Module 3
if (( $mod3 == 1 ))
then


t="$(date)" 
echo "...starting mapping at $t..."

conda activate env_MAPmiR

# BOWTIE OPTIONS
#
# -n 1 = default alignment mode and allow 1 mis-matches within that region 
# -l = defines length of "seed region"
# --norc = don't attempt to align against reverse complement
# --best; --strata = report only best alignement
# -m 1 = suppress reads mapping to more than 1 region
# --threads = number of parallel cores being used
# --un = filename for all unaligned reads
# --sam = SAM output
#
# NB 2> directs standard error file to output (appropriate for log)

# NB using -n 1 results in significantly more matches than -n 0, which is usally preferred for miRNA mapping.  Potentially justifiable given different breeds of Dog?  
# can omit genome mapping to save time


for f in trim/*.fq; do

(bowtie -n 1 -l 32 --norc --best --strata -m 1 --threads $set_cores \
  -x $set_genomes_dir/dog \
  $f \
  --un ${f}.dog.unaligned.fq \
  --sam ${f}.dog.sam) \
  2>${f}.dog.log

(bowtie -n 1 -l 32 --norc --best --strata -m 1 --threads $set_cores \
  -x $set_genomes_dir/mouse \
  ${f}.dog.unaligned.fq \
  --un ${f}.mouse.unaligned.fq \
  --sam ${f}.mouse.sam) \
  2>${f}.mouse.log
  
(bowtie -n 1 -l 32 --norc --best --strata -m 1 --threads $set_cores \
  -x $set_genomes_dir/human \
  ${f}.mouse.unaligned.fq \
  --un ${f}.human.unaligned.fq \
  --sam ${f}.human.sam) \
  2>${f}.human.log

done

mkdir map; mv trim/*dog* trim/*mouse* trim/*human* map


# sort > index > deduplicate > index > count

parallel "samtools sort {} > {.}.bam" ::: map/*.sam 
parallel "samtools index {}" ::: map/*.bam 


# ideally use --method=directional but memory usage excessive

parallel "umi_tools dedup \
 --method=unique \
 -I {} \
 -S {.}.deduplicated.bam \
 -L {.}.deduplicated.log" ::: map/*.bam
 
mkdir dedup; mv map/*.dedup* dedup

#### idxstats
#
# column 1 = gene name
# column 3 = read count
# head -n -1 to lose final empty line

parallel "samtools index {}" ::: dedup/*.deduplicated.bam
parallel "samtools idxstats {} | \
          cut -f 1,3 | \
          head -n -1 \
          > {.}.txt" ::: dedup/*.deduplicated.bam

mkdir counts; mv dedup/*.txt counts
mv trim/lib_size.txt counts

conda deactivate

fi


#### Close
t="$(date)"
echo "...finished at $t."
echo ""
echo ""
echo "Output files:"
echo "----------------"
ls -ho counts/*.txt
