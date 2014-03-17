#!/usr/bin/python
# -*- coding: utf-8 -*-
# Fonctions du tp2 du cours IMN530 - H14 - Jérémie Fouquet et Vincent Méthot

"""Doc string tp2"""

import numpy as np
import nibabel as nib
import scipy.interpolate
import Image
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors as mpc
import warnings
import time

# Enlever à la fin
import pdb

def JointHist(I, J, nbin=256, normIm=False):
    """Calcule l'histogramme conjoint de deux images de même taille (I et J)
    en divisant leur intervalle de valeurs en 'nbin' sous-intervalles

    Paramètres
    ----------
    I et J: Images (2D) en format NifTi-1, jpg ou png.
    nbin: int, optionnel. Le nombre de bins pour le calcul de
    l'histogramme. 256 par défaut.
    normIm: logical, optionnel. Si True, les images sont normalisées avant le
    calcul de l'histo. Ainsi, la plus basse (haute) valeur de chaque image
    correspond au preimer (dernier) bin. Si mis à True, peut créer des
    artéfacts dus aux arrondissements. Défaut: False.

    Exemple
    -------

    >>> H = tp2.JointHist('../Data/I4.jpg', '../Data/J4.jpg')"""
    # Ouverture et mises en forme des images en fonction des options entrées et
    # du format des images
    I = openImage(I)
    J = openImage(J)
    areImInt = (np.issubdtype(I.dtype, np.integer) &
                np.issubdtype(J.dtype, np.integer))
    if areImInt:
        # Tests sur le nombre de bins et le nombre de valeurs possibles dans
        # les images.
        maxNBytes = np.maximum(I.dtype.itemsize, J.dtype.itemsize)
        maxNValue = 2 ** (8 * maxNBytes)
        if nbin < maxNValue and not(normIm):
            warnings.warn("nbin < le nombre de valeurs possibles dans au moins "
                          "une des deux  images. Pour tracer l'histogramme, "
                          "nous devons normaliser les images, ce qui peut "
                          "introduire des artéfacts.")
            normIm = True
        minNBytes = np.minimum(I.dtype.itemsize, J.dtype.itemsize)
        minNValue = 2 ** (8 * minNBytes)
        if nbin > minNValue and not(normIm):
            warnings.warn("nbin > nombre de valeurs possibles dans au moins "
                          "une image et aucune normalisation demandée. "
                          "Nous allons normaliser les images, car sinon les "
                          "bins extrêmes ne correspondront à aucune valeur "
                          "possibles.")
            normIm = True
    elif not(normIm):
        warnings.warn("Les images ne contiennent pas des entiers. Malgré "
                      "l'option normIm=False, nous allons donc les normaliser.")
        normIm = True

#
    if normIm:
        I = np.round(normalizeIm(I) * (nbin - 1)).astype(int)
        J = np.round(normalizeIm(J) * (nbin - 1)).astype(int)

    # Impression d'info et déclaration de l'histogramme
    print I.max(), I.min(), I.std()
    print J.max(), J.min(), J.std()
    H = np.zeros([nbin, nbin], dtype=int)
    print H.shape


    # À faire: s'assurer que les deux images aient les même dimensions
    # (interpolation)
    # i, j = np.meshgrid(range(nbin), range(nbin))

    # Calcul de l'histogramme
    # Pas efficace
    # for i in range(nbin):
    #     for j in range(nbin):
    #         H[i, j] = ((I == i) & (J == j)).sum()

    # Ne va fonctionner que pour les images 2D
    for x in range(I.shape[0]):
        for y in range(I.shape[1]):
            H[I[x, y], J[x, y]] += 1

    return H


def verifSommeHisto(I, J, nbin=256):
    """Vérifie que l'histogramme conjoint de 2 images de même taille comprend
    bien, lorsque sommé sur toutes ses bins, le nombre de pixels d'une image
    Retourne True si la condition est vérifiée, False sinon.

    Paramètresimport griddata
    ----------
    I et J: Images (2D) en format NifTi-1, jpg ou png.
    nbin: int, optionnel. Le nombre de bins pour le calcul de l'histogramme.
256 par défaut.

    Exemple
    -------
    ­>>> tp2.verifSommeHisto('../Data/I1.png','../Data/J1.png')
    """

    jointHist = JointHist(I, J, nbin)
    I = openImage(I)
    return I.size == jointHist.sum()


def pltJointHist(I, J, nbin=256, normIm=False, colorMap = 'jet', cLimFrac = [0,1], useLogNorm = False):
    """Affiche l'histogramme conjoint des images I et J avec l'origine située
    dans le coin inférieur gauche.

    Paramètres
    ----------
    I, J, nbin et normIm: paramètres passés directement, dans cet ordre, à
    JointHist. Voir l'aide de JointHist pour plus d'info.
    colorMap: string, optionnel. Colormap de l'histogramme affiché.
    cLimFrac: liste 2x1, optionnelle. Spécifie les limites pour la mise à
l'echelle de l'image dans la colormap. Se specifie en terme de fraction du min
et du max de l'image. Défaut: [0,1].
    useLogNorm: logical, optionnel. Si True, une échelle log est utilisée pour
l'intensité de l'image. Défaut: False.

    Exemple
    -------
    ­>>> tp2.pltJointHist('../Data/I3.jpg','../Data/J3.jpg',cLimFrac = [0,0.0005])
    ­>>> tp2.pltJointHist('../Data/I3.jpg','../Data/J3.jpg',useLogNorm=True)
    """

    jointHist = JointHist(I, J, nbin, normIm)
    mainFig = plt.figure('IMN530 - Histo conjoint')
    mainFig.clf()
    minMax = [jointHist.min(), jointHist.max()]
    customClim = [a * b for a, b in zip(cLimFrac, minMax)]
    # Préparation de la normalisation logarithmique si demandé
    if not useLogNorm:
        customNorm = None
    else:
        if customClim[0] == 0:
            customClim[0] = 1e-15
        customNorm = mpc.LogNorm(vmin=customClim[0],vmax=customClim[1],clip=True)
    # Affichage de l'image
    imAxes = plt.imshow(jointHist, cmap=colorMap, clim = customClim,
        norm=customNorm, interpolation="none")
    imAxes.get_axes().invert_yaxis()
    #imAxes.get_axes().set_xlabel("Intensites de J")
    #imAxes.get_axes().set_ylabel("Intensites de I")
    plt.colorbar()
    plt.draw()
    plt.show(block=False)

    return mainFig, imAxes


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

    return CR

def IM(I, J, nbin=256):
    """Calcule l'information mutuelle entre 2 images de même taille

    Exemple
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
    IM = (H.astype(float) / H.sum() * log(H.sum() * H.astype(float) /
        (Hi * Hj * Hk))).sum()

    print "Information mutuelle:", IM
    return IM

def trans_rigide(theta, omega, phi, p, q, r):
    """Renvoie la matrice de transformation rigide en coordonnées homogènes
    -------------------------------------------------------------------
    theta:         angle de rotation autour de x
    omega:         angle de rotation autour de y
    phi:         angle de rotation autour de z
    (p, q, r):     vecteur de translation

    Voir aussi
    ----------
    gille_test: test d'une transformation à l'aide d'une grille de points"""

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
    s:             rapport de l'homothétie
    theta:         angle de rotation autour de x
    omega:         angle de rotation autour de y
    phi:         angle de rotation autour de z
    (p, q, r):     vecteur de translation

    Voir aussi
    ----------
    gille_test: test d'une transformation à l'aide d'une grille de points"""

    Transformation = np.matrix(np.diag([s, s, s, 1])) * trans_rigide(theta, omega, phi, p, q, r)
    return Transformation

def translation(I, p, q, interpMeth="linear"):
    """Retourne une nouvelle image correspondant à la translatée
de l'image 'I' par le vecteur t = (p, q) (p et q doivent être des float)

La gestion de l'interpolation est effectuée par scipy.interpolate

Exemple
-------

>>> J = tp2.translation('../Data/I1.png', 4.5, 6.7)"""

    I = openImage(I)
    xi = np.mgrid[:I.shape[0], :I.shape[1]]
    points = xi.astype(float)
    points[0,...] += p
    points[1,...] += q

    J = scipy.interpolate.griddata((points[0].ravel(), points[1].ravel()),
        I.ravel(), (xi[0], xi[1]), method=interpMeth, fill_value=0)

    return J


def rec2dtrans(I, J, stepSize=1e-7, pqConstCptMax=10, minDeltaPq=0.01,
               nItMax=10000, showEvo=False):
    """Recalage 2D minimisant la SSD grâce à une descente de gradient à pas
    fixe en considérant uniquement les translations. L'énergie SSD
    correspondant à chaque état est sauvegardée. L'image I est translatée pour
    correspondre à l'image J.

    Paramètres
    ----------
    I, J: Images (2D) en format NifTi-1, jpg, png, ou ndarray.
    stepSize: float, optionnel. Pas de gradient constant.
    smallGradCptMax: int, optionnel. Si la norme du gradient < constante
    déterminée à partir de stepSize pendant smallGradCptMax itérations, la
    descente s'arrete.'
    nItMax: int, optionnel. Nombre maximal d'itérations effectuées
    showEvo: logical, optionnel. Si True, l'évolution de l'image modifiée est
    montrée sur une figure.

    Retour
    ------
    ITrans: Image I translatée recalée (du moins en théorie) avec J
    allSsd: ndarray,

    Exemple
    -------
    ­>>> tp2.rec2dtrans("../Data/BrainMRI_1.jpg", "../Data/BrainMRI_2.jpg",
        showEvo=True)
    """

    # Descente de gradient à pas fixe.
    I = openImage(I)
    J = openImage(J)
    actualPq = np.array([0, 0])
    grad = np.array([np.inf, np.inf])
    # Constantes pour la condition d'arret.
    pqConstCpt = 0
    prevPq = np.array([0, 0])

    # (si norm(grad) < smallestGrad,pour
    # smallGradCptMax fois de suite, la fonction est arretée. Si nIt sont
    # faites, la fonction est aussi arretée, mais avec un warning.)
    #smallestGrad = stepSize * 1e11  # Déterminé empiriquement pour l'instant
    #smallGradCpt = 0  # Initialisation du compteur

    # Array pour stocker les SSD
    allSsd = np.zeros(nItMax)

    for iIt in range(0, nItMax):
        # Calcul de l'image translatée
        ITrans = translation(I, prevPq[0], prevPq[1])
        # Calcul du gradient selon p et q (vérifier la concordance avec les
        # "x" et "y" de la fonction translation)
        IDiff = (ITrans - J)

        gradDim0, gradDim1 = np.gradient(ITrans)
        grad[0] = -2 * (IDiff * gradDim0).sum()
        grad[1] = -2 * (IDiff * gradDim1).sum()

        # Mise à jour du vecteur de translation.
        actualPq = prevPq - stepSize * grad
        print "{}{}".format("grad = ", grad)
        print "{}{}".format("[p,q] = ", prevPq)
        # Calcul et stockage de la SSD
        allSsd[iIt] = SSD(ITrans, J)
        print "{}{}".format("SSD actuel: ", allSsd[iIt])
        # Montrer l'évolution si demandé
        if showEvo:
            pltRecalage2D(ITrans, J)
            time.sleep(0.05)
        # Condition d'arret.
        # Incrementation d'un compteur si le vecteur translation ne change pas
        # plus qu'à une précision de minDeltaPq pixel.
        if np.linalg.norm(actualPq - prevPq) < minDeltaPq:
            pqConstCpt += 1
        else:
            pqConstCpt = 0
        if pqConstCpt > pqConstCptMax:
            convReached = True
            break
        prevPq = actualPq
        # Incrémentation d'un compteur si le gradient est petit
        #if np.linalg.norm(grad) < smallestGrad:
            #smallGradCpt = smallGradCpt + 1
        #else:
            #smallGradCpt = 0
        #if smallGradCpt > smallGradCptMax:
            #convReached = True
            #break

    # Afficher un warning si la fonction s'est arretée en raison du nombre d'it
    # max
    if not(convReached):
        warnings.warn("Le recalage s'est arreté car le nombre d'itérations max"
                      " a été atteint. Le recalage peut ne pas etre optimal.")
    # Élmination des 0 à la fin de allSsd
    allSsd = np.delete(allSsd, range(np.minimum(iIt + 1, nItMax), nItMax))

    return actualPq, ITrans, allSsd


def rec2drot(I, J, stepSize=1e-7, aConstCptMax=10, minDeltaA=0.001,
             nItMax=10000, showEvo=False):
    """Recalage 2D minimisant la SSD grâce à une descente de gradient à pas
    fixe en considérant uniquement les rotations. L'énergie SSD
    correspondant à chaque état est sauvegardée. L'image I est translatée pour
    correspondre à l'image J.

    Paramètres
    ----------
    I, J: Images (2D) en format NifTi-1, jpg, png, ou ndarray.
    stepSize: float, optionnel. Pas de gradient constant.
    smallGradCptMax: int, optionnel. Si la norme du gradient < constante
    déterminée à partir de stepSize pendant smallGradCptMax itérations, la
    descente s'arrete.'
    nItMax: int, optionnel. Nombre maximal d'itérations effectuées
    showEvo: logical, optionnel. Si True, l'évolution de l'image modifiée est
    montrée sur une figure.

    Retour
    ------
    ITrans: Image I rotatée recalée (du moins en théorie) avec J
    allSsd: ndarray,

    Exemple
    -------
    ­>>> tp2.rec2drot("../Data/BrainMRI_1.jpg", "../Data/BrainMRI_3.jpg",
        showEvo=True)
    """

    # Descente de gradient à pas fixe.
    I = openImage(I)
    J = openImage(J)
    prevA = 0
    grad = np.inf
    # Constantes pour la condition d'arret.
    aConstCpt = 0
    # (si norm(grad) < smallestGrad,pour
    # smallGradCptMax fois de suite, la fonction est arretée. Si nIt sont
    # faites, la fonction est aussi arretée, mais avec un warning.)
    #smallestGrad = stepSize * 1e11  # Déterminé empiriquement pour l'instant
    #smallGradCpt = 0  # Initialisation du compteur



    # Array pour stocker les SSD
    allSsd = np.zeros(nItMax)

    for iIt in range(0, nItMax):
        # Calcul de l'image translatée
        IRot = rotation(I, prevA)
        # Calcul du gradient selon a
        IDiff = (IRot - J)
        coord = np.mgrid[:J.shape[0], :J.shape[1]]
        gradDim0, gradDim1 = np.gradient(ITrans)
        IDeriv = gradDim0 * (-np.sin(a) * coord[0, ...]
                             - np.cos(a) * coord[1, ...]) \
                 + gradDim1 * (np.cos(a) * coord[0, ...]
                               - np.sin(a) * coord[1, ...])
        grad = 2 * (IDiff * IDeriv).sum()

        # Mise à jour de l'angle de rotation.
        actualA = prevA - stepSize * grad
        print "{}{}".format("grad = ", grad)
        print "{}{}".format("[p,q] = ", actualA)
        # Calcul et stockage de la SSD
        allSsd[iIt] = SSD(ITrans, J)
        print "{}{}".format("SSD actuel: ", allSsd[iIt])
        # Montrer l'évolution si demandé
        if showEvo:
            pltRecalage2D(IRot, J)
            time.sleep(0.05)
        # Condition d'arret.
        # Incrementation d'un compteur si l'angle ne change pas
        # plus qu'à une précision de minDeltaA pixel.
        if (actualA - prevA) < minDeltaA:
            aConstCpt += 1
        else:
            aConstCpt = 0
        if aConstCpt > pqConstCptMax:
            convReached = True
            break
        prevA = actualA

    # Afficher un warning si la fonction s'est arretée en raison du nombre d'it
    # max
    if not(convReached):
        warnings.warn("Le recalage s'est arreté car le nombre d'itérations max"
                      " a été atteint. Le recalage peut ne pas etre optimal.")
    # Élmination des 0 à la fin de allSsd
    allSsd = np.delete(allSsd, range(np.minimum(iIt + 1, nItMax),
                                                nItMax))

    return pq, ITrans, allSsd


def rotation(I, theta):
    """Application d'une rotation d'angle 'theta' et de centre (0, 0)
    (coin supérieur gauche) à l'image 'I'"""

def rec2dpasfixe(I, J):
    """Recalage 2D minimisant la SSD par une descente de gradient.
    Considère l'ensemble des transformations rigides."""

def rec2doptimize(I, J):
    """Recalage 2D minimisant la SSD par une descente de gradient optimisée.
    Considère l'ensemble des transformations rigides."""



def openImage(I):
    """Ouvre des images au format jpeg, png et NifTI et les retourne en numpy
    array. (si complexe, prend la valeur absolue)"""

    if isinstance(I, str):
        if (I[-7:] == '.nii.gz') | (I[-4:] == '.nii'):
            #J = np.asarray(nib.load(I).get_data(), dtype=float)
            J = np.asarray(nib.load(I).get_data())
        elif (I[-4:] == '.jpg') | (I[-5:] == '.jpeg') | (I[-4:] == '.png'):
            #J = np.asarray(Image.open(I), dtype=float)
            J = np.asarray(Image.open(I))
        else:
            print "Formats d'image acceptés: .nii, .nii.gz, .jpg, .png, .jpeg"
    else:
        #J = np.abs(np.asarray(I)).astype(float)
        J = np.abs(np.asarray(I))

    # Normalisation de l'image
    #if rescaleIm:
    #    J = (J - J.min()) / (J - J.min()).max()

    return J


def normalizeIm(I):
    """ Transforme une image en float et rescale ses valeurs entre 0 et 1.'"""
    I = I.astype(float)
    return (I - I.min()) / (I - I.min()).max()


def grille_test(transformation, xmax = 6, ymax = 6, zmax = 2):
    """Test des transformations à l'aide d'une grille de points.
---------------------------------------------------------
transformation: matrice de transformation 3D en coordonnées homogènes
xmax, ymax, zmax: limites de la grille de point (min à 0)

Exemple:
--------
>>> tp2.grille_test([[2, 0, 0, 5.5], [0, 3, 0, 4.5], [0, 0, 1, 4.5],
[0, 0, 0, 1]])
>>> tp2.grille_test(tp2.M3)"""

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xis, yis, zis = np.mgrid[:xmax,:ymax,:zmax]
    pis = np.matrix([xis.ravel(), yis.ravel(), zis.ravel(),
        ones(xis.shape).ravel()]).T
    ax.scatter(xis, yis, zis, c='b')

    pfs = pis * transformation
    xfs, yfs, zfs, uns = array(pfs.T)
    ax.scatter(xfs, yfs, zfs, c='r')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.axis('tight')

    plt.show()


def pltRecalage2D(I, J):
    """Affiche 2 images de même taille ainsi que la différence entre ces deux
    images"""

    I = openImage(I)
    J = openImage(J)

    mainFig = plt.figure("IMN530 - Recalage images 2D")
    currAx = plt.subplot(131)
    plt.imshow(I)
    plt.title("Image 1")
    currAx.xaxis.set_visible(False)
    currAx.yaxis.set_visible(False)

    currAx = plt.subplot(132)
    plt.imshow(J)
    plt.title("Image 2")
    currAx.xaxis.set_visible(False)
    currAx.yaxis.set_visible(False)

    currAx = plt.subplot(133)
    plt.imshow(np.abs(I - J))
    plt.title("Diff abs")
    currAx.xaxis.set_visible(False)
    currAx.yaxis.set_visible(False)

    plt.draw()
    plt.show(block=False)

    return mainFig


def num3b(stepSize=1e-7, smallGradCptMax=10, nItMax=10000, showEvo=True):
    """Fonction qui teste rec2dtrans pour 3 translations aléatoires différentes
    pour répondre à la question 3.b. du tp.

    Paramètres
    ----------
    stepSize, smallGradCptMax, nItMax et showEvo: Même param que pour
    rec2dtrans. Voir aide de cette fonction.

    Exemple
    -------
    regPq, initPq, allSsd = tp2.num3b()
    """

    I = openImage("../Data/BrainMRI_1.jpg")
    # Array to store the 3 initial and post-registration translations
    initPq = np.zeros([2, 3])
    regPq = np.zeros([2, 3])
    # Make the whole process reproducible
    np.random.seed(1)
    # Make array to store ssd
    allSsd = [0, 0, 0]
    # Do registration for the 3 translations
    for iPq in range(0, 3):
        initPq[:, iPq] = (np.random.rand(2) * np.array(I.shape) - 0.5) / 4
        ITrans = translation(I, initPq[0, iPq], initPq[1, iPq])
        regPq[:, iPq], dum, allSsd[iPq] = rec2dtrans(ITrans, I, stepSize=1e-7,
                                          nItMax=10000, showEvo=True)

    return regPq, initPq, allSsd



# Transformations de la question 3d

M1 = np.matrix([[0.9045, -0.3847, -0.1840, 10.0000],
                [0.2939,  0.8750, -0.3847, 10.0000],
                [0.3090,  0.2939,  0.9045, 10.0000],
                [0.0000,  0.0000,  0.0000,  1.0000]])

M2 = np.matrix([[0.0000, -0.2598,  0.1500, -3.0000],
                [0.0000, -0.1500, -0.2598,  1.5000],
                [0.3000,  0.0000,  0.0000,  0.0000],
                [0.0000,  0.0000,  0.0000,  1.0000]])

M3 = np.matrix([[ 0.7182, -1.3727, -0.5660,  1.8115],
                [-1.9236, -4.6556, -2.5512,  0.2873],
                [-0.6426, -1.7985, -1.6285,  0.7404],
                [ 0.0000,  0.0000,  0.0000,  1.0000]])


