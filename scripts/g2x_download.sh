#!/bin/bash
# G2-X external cohort: 10x Visium "Parent_Visium_Human_BreastCancer" (v1.2.0), public no-auth.
# URLs live-verified HTTP 200 on 2026-09-27 and again 2026-09-28 (content-lengths below).
# Lost once to a sandbox wipe because they were not committed; do not let that happen again.
set -e
mkdir -p /home/sandbox/g2x && cd /home/sandbox/g2x
BASE='https://cf.10xgenomics.com/samples/spatial-exp/1.2.0/Parent_Visium_Human_BreastCancer'
curl -O "$BASE/Parent_Visium_Human_BreastCancer_filtered_feature_bc_matrix.h5"   # 21,452,440 bytes
curl -O "$BASE/Parent_Visium_Human_BreastCancer_spatial.tar.gz"                  # 11,759,586 bytes
curl -o image.tif "$BASE/Parent_Visium_Human_BreastCancer_image.tif"             # 1,763,006,484 bytes
tar xzf Parent_Visium_Human_BreastCancer_spatial.tar.gz
