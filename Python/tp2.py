#!/usr/bin/python
# -*- coding: utf-8 -*-
# Fonctions du tp2 du cours IMN530 - H14 - Jérémie Fouquet et Vincent Méthot

import numpy as np
import nibabel as nib
import Image
from pylab import *

def JointHist(I, J, nbin=256):
	"""Calcule l'histogramme conjoint de deux images de même taille (I et J)
	en divisant leur intervalle de valeurs en 'nbin' sous-intervalles

	Example
	-------

	>>> H = tp2.JointHist('../Data/I4.jpg', '../Data/J4.jpg')"""

	I = ((nbin-1) * openImage(I)).astype(int)
	print I.max(), I.min(), I.std()
	J = ((nbin-1) * openImage(J)).astype(int)
	print J.max(), J.min(), J.std()
	H = np.zeros([nbin, nbin], dtype=int)
	print H.shape

	# À faire: s'assurer que les deux images aient les même dimensions (interpolation)

	# i, j = np.meshgrid(range(nbin), range(nbin))

	# Pas efficace
	# for i in range(nbin):
	# 	for j in range(nbin):
	# 		H[i, j] = ((I == i) & (J == j)).sum()

	# Ne va fonctionner que pour les images 2D

	for x in range(I.shape[0]):
		for y in range(I.shape[1]):
			H[I[x,y], J[x,y]] += 1

	return H

def SSD(I, J, nbin=256):
	"""Calcule la somme des différences au carré entre 2 images (I et J)
	de même taille"""

	# À partir de l'histogramme
	H = JointHist(I, J, nbin)
	i, j = np.meshgrid(range(nbin), range(nbin))
	SSD = (H*(i - j)**2).sum()

	# Sans histogramme

	return SSD

def CR(I, J, nbins=256):
	"""Calcule le coefficient de corrélation entre 2 images (I et J)
	de même taille"""

	# À partir de l'histogramme
	H = JointHist(I, J, nbin)
	i, j = np.meshgrid(range(nbin), range(nbin))
	meanI = 1. / H.sum() * (H * i).sum()
	meanJ = 1. / H.sum() * (H * j).sum()
	covariance = (H * (i - meanI) * (j - meanJ)).sum()
	autocovI = 1. / H.sum() * (H * i**2 - meanI**2).sum()
	autocovJ = 1. / H.sum() * (H * j**2 - meanJ**2).sum()
	CR = covariance / np.sqrt(autocovI * autocovJ)

	# Sans histogramme

	return CR

def IM(I, J, nbin=256):
	"""Calcule l'information mutuelle entre 2 images de même taille

	Example
	-------

	>>> IM = tp2.JointHist('../Data/I4.jpg', '../Data/J4.jpg')"""

	# À partir de l'histogramme
	H = JointHist(I, J, nbin)

	# Broadcasting!
	Hi = np.empty(list(H.shape))
	Hj = np.empty(list(H.shape))
	Hi[:] = H.sum(0)
	Hj[:] = H.sum(1)

	if H.ndim == 3:
		Hk = np.empty(list(H.shape)) 
		Hk[:] = H.sum(2)
		Hk[Hk == 0] = Hk.sum()
	else:
		Hk = 1

	# On remplace les 0 par des 1 dans l'histogramme
	# Il reste plein de choses à optimiser...

	H[H == 0] = H.sum()
	Hi[Hi == 0] = Hi.sum()
	Hj[Hj == 0] = Hj.sum()
	IM = (H.astype(float) / H.sum() * log(H.sum() * H.astype(float) / (Hi * Hj * Hk))).sum()

	print "Information mutuelle:", IM	
	return IM

def trans_rigide(theta, omega, phi, p, q, r):
	"""Renvoie la matrice de transformation rigide en coordonnées homogènes
	-------------------------------------------------------------------
	theta: 		angle de rotation autour de x
	omega: 		angle de rotation autour de y
	phi: 		angle de rotation autour de z
	(p, q, r): 	vecteur de translation"""

	T = np.matrix([[1, 0, 0, p], [0, 1, 0, q], [0, 0, 1, r], [0, 0, 0, 1]])
	Rx = np.matrix([[1, 0, 0, 0], [0, np.cos(theta), -np.sin(theta), 0], 
		[0, np.sin(theta), np.cos(theta), 0], [0, 0, 0, 1]])
	Ry = np.matrix([[np.cos(omega), 0, -np.sin(omega), 0], [0, 1, 0, 0], 
		[np.sin(omega), 0, np.cos(omega), 0], [0, 0, 0, 1]])
	Rz = np.matrix([[np.cos(phi), -np.sin(phi), 0, 0], [np.sin(phi), np.cos(phi), 0, 0], 
		[0, 0, 1, 0], [0, 0, 0, 1]])
	Transformation = T * Rz * Ry * Rx

	return Transformation

def similitude(s, theta, omega, phi, p, q, r):
	"""Revoie la matrice de transformation rigide (+homothétie)
	en coordonnées homogènes
	-------------------------------------------------------------------
	s: 			rapport de l'homothétie
	theta: 		angle de rotation autour de x
	omega: 		angle de rotation autour de y
	phi: 		angle de rotation autour de z
	(p, q, r): 	vecteur de translation"""

def translation(I, p, q):
	"""Retourne une nouvelle image correspondant à la translatée
	de l'image 'I' par le vecteur t = (p, q) (p et q doivent être des float)

	La gestion de l'interpolation est effectuée par [module d'interpolation]"""

def rec2dtrans(I, J):
	"""Recalage 2D minimisant la SSD et considérant uniquement les translations.
	L'énergie SSD correspondant à chaque état est sauvegardée."""

def rotation(I, theta):
	"""Application d'une rotation d'angle 'theta' et de centre (0, 0)
	(coin supérieur gauche) à l'image 'I'"""

def rec2drot(I, J):
	"""Recalage 2D minimisant la SSD et considérant uniquement les rotations.
	L'énergie SSD correspondant à chaque état est sauvegardée."""

def rec2dpasfixe(I, J):
	"""Recalage 2D minimisant la SSD par une descente de gradient.
	Considère l'ensemble des transformations rigides."""

def rec2doptimize(I, J):
	"""Recalage 2D minimisant la SSD par une descente de gradient optimisée.
	Considère l'ensemble des transformations rigides."""



def openImage(I):
	"""Ouvre des images au format jpeg, png et NifTI et les retourne en numpy array.
	Normalise et transforme en float array les autres types d'entrée (si complexe,
	prend la valeur absolue)"""

	if isinstance(I, str):
		if (I[-7:] == '.nii.gz') | (I[-4:] == '.nii'):
			J = np.asarray(nib.load(I).get_data(), dtype=float)
		elif (I[-4:] == '.jpg') | (I[-5:] == '.jpeg') | (I[-4:] == '.png'):
			J = np.asarray(Image.open(I), dtype=float)
		else:
			print "Formats d'image acceptés: .nii, .nii.gz, .jpg, .png, .jpeg"
	else:
		J = np.abs(np.asarray(I)).astype(float)

	# Normalisation de l'image
	J = (J - J.min()) / (J - J.min()).max()

	return J
		

