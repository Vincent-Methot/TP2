#!/usr/bin/python
# -*- coding: utf-8 -*-
# Fonctions du tp2 du cours IMN530 - H14 - Jérémie Fouquet et Vincent Méthot

import numpy as np

def JointHist(I, J, bin):
	"""Calcule l'histogramme conjoint de deux images de même taille (I et J)
	en divisant leur intervalle de valeurs en 'bin' sous-intervalles"""