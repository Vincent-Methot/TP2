#!/usr/bin/python
# -*- coding: utf-8 -*-
# Tests reliés au tp2 du cours IMN530 - H14 - Jérémie Fouquet et Vincent Méthot

import tp2
import numpy as np


# Question 3

grille = np.mgrid[:20,:20,:5]
print grille.shape


M1 = np.matrix([[0.9045, -0.3847, -0.184, 10.], [0.2939, 0.8750, -0.3847, 10.], 
	[0.3090, 0.2939, 0.9045, 10.], [0, 0, 0, 1.]])
M2 = np.matrix([[0, -0.2598, 0.1500, -3.], [0, -0.15, -0.2598, 1.5], 
	[0.3, 0, 0, 0], [0, 0, 0, 1]])
M3 = np.matrix([[0.7182, -1.3727, -0.5660, 1.8115], [-1.9236, -4.6556, -2.5512, 0.2873], 
	[-0.6426, -1.7985, -1.6285, 0.7404], [0, 0, 0, 1.]])