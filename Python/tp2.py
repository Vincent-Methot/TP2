#!/usr/bin/python
# -*- coding: utf-8 -*-
# Fonctions du tp2 du cours IMN530 - H14 - Jérémie Fouquet et Vincent Méthot

import numpy as np
import nibabel as nib
import Image

def JointHist(I, J, bin):
	"""Calcule l'histogramme conjoint de deux images de même taille (I et J)
	en divisant leur intervalle de valeurs en 'bin' sous-intervalles"""

	I = openImage(I)
	J = openImage(J)


def SSD(I, J):
	"""Calcule la somme des différences au carré entre 2 images (I et J)
	de même taille"""

def CR(I, J):
	"""Calcule le coefficient de corrélation entre 2 images (I et J)
	de même taille"""

def IM(I, J):
	"""Calcule l'information mutuelle entre 2 images de même taille"""

def trans_rigide(theta, omega, phi, p, q, r):
	"""Renvoie la matrice de transformation rigide en coordonnées homogènes
	-------------------------------------------------------------------
	theta: 		angle de rotation autour de x
	omega: 		angle de rotation autour de y
	phi: 		angle de rotation autour de z
	(p, q, r): 	vecteur de translation"""

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
	"""Ouvre des images au format jpeg, png et NifTI et les retourne en numpy asarray"""

	if I[-7:] == '.nii.gz' | I[-4:] == '.nii':
		J = np.asarray(nib.load(I).get_data())
	elif I[-4:] == '.jpg' | I[-4:] == '.png' | I[-5] == '.jpeg':
		J = np.asarray(Image.open(I))
	else:
		print "Formats d'image acceptés: .nii, .nii.gz, .jpg, .png, .jpeg"

	return J
		
# Test par J. Fouquet
