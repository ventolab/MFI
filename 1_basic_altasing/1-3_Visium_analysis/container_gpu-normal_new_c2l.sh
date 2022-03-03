#!/usr/bin/env bash

bsub -q gpu-normal -M200000 \
  -G team292 \
  -R"select[mem>200000] rusage[mem=200000, ngpus_physical=1.00] span[hosts=1]"  \
  -gpu "mode=shared:j_exclusive=yes" -Is \
  /software/singularity-v3.5.3/bin/singularity exec \
  --no-home  \
  --nv \
 -B /nfs,/lustre \
  /nfs/cellgeni/singularity/images/cell2location-4b1c05c.sif \
  /bin/bash -c "HOME=$(mktemp -d) jupyter notebook --notebook-dir=/nfs/team292/aa22/spatial_deconvolution_related/202008_real_MFI_Visium/202111_new_analysis_all_samples --NotebookApp.token='cell2loc' --ip=0.0.0.0 --port=1234 --no-browser --allow-root"
