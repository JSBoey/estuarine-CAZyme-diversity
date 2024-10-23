#!/bin/bash -e
#SBATCH --job-name=drep_98ani
#SBATCH --account=ga02676
#SBATCH --time=3:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=24
#SBATCH --output=slurm_out/%j.%x.out
#SBATCH --error=slurm_err/%j.%x.err
#SBATCH --partition=milan

module purge
module load drep/3.4.2-gimkl-2022a-Python-3.10.5

genome_dir=data/00-mags-unrefined-rename-out
work_dir=data/drep_98ani_out
checkm_csv=data/checkm_for_drep.csv

dRep dereplicate \
  --processors ${SLURM_CPUS_PER_TASK} \
  --genomes ${genome_dir}/*.fa \
  --genomeInfo ${checkm_csv} \
  --S_algorithm gANI \
  --completeness 0 \
  --contamination 100 \
  --P_ani 0.90 \
  --S_ani 0.98 \
  ${work_dir}
