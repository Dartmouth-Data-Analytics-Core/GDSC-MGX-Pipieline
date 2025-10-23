#!/usr/bin/env bash

#SBATCH --job-name=downsample
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=downsample_%j.out

#----- START
echo -e "#-----------------------------------------------------------#"
echo "Starting job: $SLURM_JOB_NAME (Job ID: $SLURM_JOB_ID)"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo -e "#-----------------------------------------------------------#"

#----- Activate conda environment
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/seqtk_env

#----- Specify input directory and output directory
INPUT="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/martinez/GDSC-MGX-Pipeline/data"
OUTPUT="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/martinez/GDSC-MGX-Pipeline/data"
cd "$INPUT"

#----- Begin the downsampling
seqtk sample -s 42 193258_01_S1_R1_001.fastq.gz 0.2 | gzip > "$OUTPUT/193258_01_S1_R1_001.fastq.gz"
seqtk sample -s 42 193258_01_S1_R2_001.fastq.gz 0.2 | gzip > "$OUTPUT/193258_01_S1_R2_001.fastq.gz"
seqtk sample -s 42 193258_02_S2_R1_001.fastq.gz 0.2 | gzip > "$OUTPUT/193258_02_S2_R1_001.fastq.gz"
seqtk sample -s 42 193258_02_S2_R2_001.fastq.gz 0.2 | gzip > "$OUTPUT/193258_02_S2_R2_001.fastq.gz"

echo "Downsampling completed."
echo "End time: $(date)"
