#! /bin/bash
# ANTS.sh

# Estimation de la transformation non-linéaire
ANTS 3 -m MI[t1.nii, flair.nii, 1, 32] -o flair_to_t1.nii.gz -i 30x20x10 -r Gauss[3,1] -t Elast[3]
# Application de la transformation non-linéaire
WarpImageMultiTransform 3 flair.nii flair_to_t1.nii.gz -R t1.nii flair_to_t1Warp.nii.gz flair_to_t1Affine.txt

# Recalage linéaire
flirt -dof 6 -cost mutualinfo -in flair.nii -out flair_lineaire.nii.gz -ref t1.nii 