#! /bin/bash
# ANTS.sh

ANTS 3 -m CC[t1.nii, flair.nii, 1, 2] -r Gauss[3,0] -t SyN[0.25] 90x100x80 -o output.nii.gz
WarpImageMultiTransform 3 flair.nii flair_to_t1.nii.gz flair_to_t1Warp.nii.gz flair_to_t1Affine.txt