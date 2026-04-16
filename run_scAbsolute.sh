snakemake --cores 16 --snakefile workflow/Snakefile_absolute \
    --workflow-profile workflow/profiles/epyc --profile slurm\
    --software-deployment-method conda apptainer --use-conda --use-singularity --rerun-incomplete \
    --singularity-args "-B /home/${USER}/.cache -B /home/${USER}/Pipelines/scAbsolute:/opt/scAbsolute -B /mnt/scratchc/fmlab/yates02/Pipelines/snakemake-illumina-alignment-scMulti/results/bam/bwamem/"





