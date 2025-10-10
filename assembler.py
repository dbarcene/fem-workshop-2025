#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : assembler.py
# Author            : David Barcene <david.barcene@utp.ac.pa>
# Date              : 26.06.2022
# Last Modified Date: 10.10.2025
# Last Modified By  : David Barcene <david.barcene@utp.ac.pa>

import numpy as np


def assemble(NL, EL, epsilon_x, epsilon_y, beta, rho):
    """
    Discretizes a rectangular area into finite elements.

    Args:
        NL: numpy array with the positions of the nodes.
        EL: numpy array with the node numbers that make up each element.
        epsilon_x: 
        epsilon_y: 
        beta: 
        rho: Global charge distribution. Numpy array for the charge distribution. 
             values: 1 charge, 0 no charge.

    Returns:
        K:
        b:
    """

    NoN = len(NL[:, 0])
    NoE = len(EL[:, 0])

    K = np.zeros([NoN, NoN])
    b = np.zeros(NoN)
    Me = np.zeros([3, 3])  # M matrix per element
    Te = np.zeros([3, 3])  # T matrix per element
    rho_e = np.zeros(3)  # rho matrix per element

    for e in range(NoE):
        x21 = NL[EL[e, 1]-1, 0] - NL[EL[e, 0]-1, 0]
        x31 = NL[EL[e, 2]-1, 0] - NL[EL[e, 0]-1, 0]
        x32 = NL[EL[e, 2]-1, 0] - NL[EL[e, 1]-1, 0]
        x13 = NL[EL[e, 0]-1, 0] - NL[EL[e, 2]-1, 0]
        y12 = NL[EL[e, 0]-1, 1] - NL[EL[e, 1]-1, 1]
        y21 = NL[EL[e, 1]-1, 1] - NL[EL[e, 0]-1, 1]
        y31 = NL[EL[e, 2]-1, 1] - NL[EL[e, 0]-1, 1]
        y23 = NL[EL[e, 1]-1, 1] - NL[EL[e, 2]-1, 1]

        # Area of constant triangular element
        Ae = 0.5*(x21*y31 - x31*y21)

        # Evaluation of de Me matrix
        Me[0, 0] = -(epsilon_x[e]*y23**2 + epsilon_y[e]*x32**2)/(4*Ae)
        Me[0, 1] = -(epsilon_x[e]*y23*y31 + epsilon_y[e]*x32*x13)/(4*Ae)
        Me[1, 0] = Me[0, 1]
        Me[0, 2] = -(epsilon_x[e]*y23*y12 + epsilon_y[e]*x32*x21)/(4*Ae)
        Me[2, 0] = Me[0, 2]
        Me[1, 1] = -(epsilon_x[e]*y31**2 + epsilon_y[e]*x13**2)/(4*Ae)
        Me[1, 2] = -(epsilon_x[e]*y31*y12 + epsilon_y[e]*x13*x21)/(4*Ae)
        Me[2, 1] = Me[1, 2]
        Me[2, 2] = -(epsilon_x[e]*y12**2 + epsilon_y[e]*x21**2)/(4*Ae)

        # Evaluation of de M matrix
        Te[0, 0] = beta[e]*Ae/6
        Te[0, 1] = beta[e]*Ae/12
        Te[1, 0] = beta[e]*Ae/12
        Te[0, 2] = beta[e]*Ae/12
        Te[2, 0] = beta[e]*Ae/12
        Te[1, 1] = beta[e]*Ae/6
        Te[1, 2] = beta[e]*Ae/12
        Te[2, 1] = beta[e]*Ae/12
        Te[2, 2] = beta[e]*Ae/6

        # Me + Te = Kij
        Ke = Me + Te

        # Evaluation of vector g
        rho_e[0] = rho[EL[e][0]-1]*Ae/3
        rho_e[1] = rho[EL[e][1]-1]*Ae/3
        rho_e[2] = rho[EL[e][2]-1]*Ae/3

        for i in range(3):
            for j in range(3):
                K[EL[e, i]-1, EL[e, j]-1] = K[EL[e, i]-1, EL[e, j]-1] + Ke[i, j]

            b[EL[e, i]-1] = b[EL[e, i]-1] + ge[i]

    return K, b
