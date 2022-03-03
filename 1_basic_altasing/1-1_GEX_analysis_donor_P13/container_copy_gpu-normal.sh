#!/usr/bin/env bash
                    
                    
bsub -q gpu-normal -M200000 \
  -G team292 \
  -R"select[mem>200000] rusage[mem=200000, ngpus_physical=1.00] span[hosts=1]"  \
  -gpu "mode=shared:j_exclusive=yes" -Is \
  /software/singularity-v3.6.4/bin/singularity exec \
  --no-home  \
  --nv \
   -B /nfs,/lustre \
   -B /lustre \
  /nfs/cellgeni/singularity/images/scvi-vento-v0.2.1.sif \
  /bin/bash -c "HOME=$(mktemp -d) jupyter notebook --notebook-dir=/nfs/team292/aa22/scVI_related/202111_upd_in_vivo_analysis --NotebookApp.token='team292' --ip=0.0.0.0 --port=1234 --no-browser --allow-root"

