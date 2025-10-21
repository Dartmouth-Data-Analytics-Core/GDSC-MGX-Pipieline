#!/bin/bash

#SBATCH --job-name=MGX
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16  
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=MGX_%j.out

#----- START
echo "#------------------------ DEBUGGING ------------------------#\n"
echo "Starting job: $SLURM_JOB_NAME (Job ID: $SLURM_JOB_ID)"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"

#----- Source conda and activate snakemake
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake

#----- Run Snakemake workflow
snakemake -s Snakefile \
	--conda-frontend conda \
	--use-conda \
	--profile cluster_profile \
	--rerun-incomplete \
	--keep-going 

#----- END
echo "End time: $(date)"