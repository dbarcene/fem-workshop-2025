#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : boundary_conditions.py
# Author            : David Barcene <david.barcene@utp.ac.pa>
# Date              : 17.07.2022
# Last Modified Date: 11.10.2025
# Last Modified By  : David Barcene <david.barcene@utp.ac.pa>
import numpy as np


def boundary_conditions(lx, ly, NoN, NL, b, K, V):
    """
    Applies Dirichlet boundary conditions.

    Args

    Returns
    """
    xn = NL[:, 0]
    yn = NL[:, 1]

    # Definition Dirichlet Boundary Conditions
    EBC_nodes = []
    EBC_vals = []
    NBC = 0  # Number of Boundary Conditions
    for i in range(NoN):
        if ((xn[i] == 0) or (xn[i] == lx) or (yn[i] == 0)):
            EBC_nodes.append(int(i))
            EBC_vals.append(V)
            NBC += 1
        elif (yn[i] == ly):
            EBC_nodes.append(int(i))
            EBC_vals.append(V)
            NBC += 1

    for i in range(NBC):
        for j in range(NoN):
            if j != EBC_nodes[i]:
                b[j] = b[j] - K[j, EBC_nodes[i]]*EBC_vals[i]

        K[:, EBC_nodes[i]] = 0
        K[EBC_nodes[i], :] = 0
        K[EBC_nodes[i], EBC_nodes[i]] = 1
        b[EBC_nodes[i]] = EBC_vals[i]
