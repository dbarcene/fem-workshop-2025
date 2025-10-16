#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : assembler.py
# Author            : David Barcene <david.barcene@utp.ac.pa>
# Date              : 26.06.2022
# Last Modified Date: 16.10.2025
# Last Modified By  : David Barcene <david.barcene@utp.ac.pa>

import numpy as np


def assemble(NL, EL, alpha_x, alpha_y, g, beta=0):
    """
    Assembles the matrix system for the general form of the Poisson ecuation
    already discretized into finite elements.

    THIS EXAMPLE ONLY USES TRIANGULAR FINITE ELEMENTS in a homogeneous pattern.
    Coarsed meshes of triangular and quadrilateral elements are not implemented.

    Args:
        NL: Numpy array with the positions of the nodes.
        EL: Numpy array with the node numbers that make up each element.
        alpha_x: Constant linear argument in the x direction for anisotropic
            medium. This argumetn multiplies the derivative in the x direction.
            Hiher order anisotropic functions are not implemented in this
            example.
        alpha_y: Constant linear argument in the y direction for anisotropic
            medium. This argumetn multiplies the derivative in the x direction.
            Hiher order anisotropic functions are not implemented in this
            example.
        beta: Constant argument or scalar function that multiplies the primary
            unknown cuantity. For simple problems this argument is zero.
        g: Numpy array with the positions of the sources.
    Returns:
        K: Assembled system matrix of dimensions [NoN, NoN], where NoN is the
            Number of Nodes in the system.
        b: Assembled system result vector of dimensions [NoN, 1], where NoN is
            the Number of Nodes in the system.
    """

    NoN = len(NL[:, 0])  # Number of Nodes
    NoE = len(EL[:, 0])  # Number of Elements

    K = np.zeros([NoN, NoN])
    b = np.zeros(NoN)
    Me = np.zeros([3, 3])  # M matrix per element
    ge = np.zeros(3)  # g matrix per element

    if beta != 0:
        Te = np.zeros([3, 3])  # T matrix per element

    for e in range(NoE):
        x21 = NL[EL[e, 1]-1, 0] - NL[EL[e, 0]-1, 0]  # x2 -x1
        x31 = NL[EL[e, 2]-1, 0] - NL[EL[e, 0]-1, 0]  # x3 -x1
        x32 = NL[EL[e, 2]-1, 0] - NL[EL[e, 1]-1, 0]  # x3 -x2
        x13 = NL[EL[e, 0]-1, 0] - NL[EL[e, 2]-1, 0]
        y12 = NL[EL[e, 0]-1, 1] - NL[EL[e, 1]-1, 1]
        y21 = NL[EL[e, 1]-1, 1] - NL[EL[e, 0]-1, 1]
        y31 = NL[EL[e, 2]-1, 1] - NL[EL[e, 0]-1, 1]
        y23 = NL[EL[e, 1]-1, 1] - NL[EL[e, 2]-1, 1]

        # Area of constant triangular element
        Ae = 0.5*(x21*y31 - x31*y21)

        # Evaluation of de Me matrix
        Me[0, 0] = -(alpha_x[e]*y23**2 + alpha_y[e]*x32**2)/(4*Ae)
        Me[0, 1] = -(alpha_x[e]*y23*y31 + alpha_y[e]*x32*x13)/(4*Ae)
        Me[1, 0] = Me[0, 1]
        Me[0, 2] = -(alpha_x[e]*y23*y12 + alpha_y[e]*x32*x21)/(4*Ae)
        Me[2, 0] = Me[0, 2]
        Me[1, 1] = -(alpha_x[e]*y31**2 + alpha_y[e]*x13**2)/(4*Ae)
        Me[1, 2] = -(alpha_x[e]*y31*y12 + alpha_y[e]*x13*x21)/(4*Ae)
        Me[2, 1] = Me[1, 2]
        Me[2, 2] = -(alpha_x[e]*y12**2 + alpha_y[e]*x21**2)/(4*Ae)

    if Te:
        # Evaluation of de M matrix
        Te[0, 0] = beta[e]*Ae/6
        Te[0, 1] = beta[e]*Ae/12
        Te[0, 2] = beta[e]*Ae/12
        Te[1, 0] = beta[e]*Ae/12
        Te[1, 1] = beta[e]*Ae/6
        Te[1, 2] = beta[e]*Ae/12
        Te[2, 0] = beta[e]*Ae/12
        Te[2, 1] = beta[e]*Ae/12
        Te[2, 2] = beta[e]*Ae/6

        # Me + Te = Kij
        Ke = Me + Te
    else:
        Ke = Me

    # Evaluation of vector g
    ge[0] = g[EL[e][0]-1]*Ae/3
    ge[1] = g[EL[e][1]-1]*Ae/3
    ge[2] = g[EL[e][2]-1]*Ae/3

    for i in range(3):
        for j in range(3):
            K[EL[e, i]-1, EL[e, j]-1] = K[EL[e, i]-1, EL[e, j]-1] + Ke[i, j]

            b[EL[e, i]-1] = b[EL[e, i]-1] + ge[i]

    return K, b
