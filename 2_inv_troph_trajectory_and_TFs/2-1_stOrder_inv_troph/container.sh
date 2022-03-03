#!/usr/bin/env bash
  
bsub -q gpu-normal -M200000 \
  -G team292 \
  -R"select[mem>200000] rusage[mem=200000, ngpus_physical=1.00] span[hosts=1]"  \
  -gpu "mode=shared:j_exclusive=yes" -Is \
  /software/singularity-v3.5.3/bin/singularity exec --nv -B /nfs,/lustre /nfs/cellgeni/singularity/images/spatialde-v0.1.sif jupyter notebook --notebook-dir=/nfs/team292/aa22/SpatialDE_colocation_model_related --port=18888 --no-browser --ip=0.0.0.0
