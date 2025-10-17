# -*- coding: utf-8 -*-
# File              : uniform_mesh.py
# Author            : David Barcene <david.barcene@utp.ac.pa>
# Date              : 10.06.2022
# Last Modified Date: 16.10.2025
# Last Modified By  : David Barcene <david.barcene@utp.ac.pa>


import numpy as np


def uniform_mesh(d1, d2, p, m, element_type='2DTR3N'):
    """
    Discretizes a rectangular area into finite elements.

    Args:
        d1: distance in the horizontal direction
        d2: distance in the vertical direction
        p: number of elements in the horizontal direction
        m: number of elements in the vertical direction
        element_type: triangular '2DTR3N' (default value), or cuadrilateral
        '2DQU4N'.

    Returns:
        NL: numpy array with the positions of the nodes.
        EL: numpy array with the node numbers that make up each element.
    """

    # Corners of the mesh
    q = np.array([[0, 0], [d1, 0], [0, d2], [d1, d2]])
    PD = 2
    NoN = (p+1)*(m+1)  # Number of Nodes

    """
    element_type: 2DQU4N
        2D: 2 dimensions
        QU: quadrilateral
        4N: 4 nodes per element
    """

    NoE = p*m          # Number of Elements
    NpE = 4            # Nodes per Element

    # Generating Nodes

    NL = np.zeros([NoN, PD])

    a = (q[1, 0] - q[0, 0])/p  # Increment in the horizontal direction
    b = (q[2, 1] - q[0, 1])/m  # Increment in the vertical direction
    n = 0  # Counter for rows in NL

    for i in range(1, m+2):
        for j in range(1, p+2):

            NL[n, 0] = q[0, 0] + (j-1)*a  # for x values
            NL[n, 1] = q[0, 1] + (i-1)*b  # for y values

            n += 1

    # Generating Elements

    EL = np.zeros([NoE, NpE])

    for i in range(1, m+1):
        for j in range(1, p+1):

            if j == 1:
                EL[(i-1)*p+j-1, 0] = (i-1)*(p+1) + j
                EL[(i-1)*p+j-1, 1] = EL[(i-1)*p+j-1, 0] + 1
                EL[(i-1)*p+j-1, 3] = EL[(i-1)*p+j-1, 0] + (p+1)
                EL[(i-1)*p+j-1, 2] = EL[(i-1)*p+j-1, 3] + 1

            else:
                # second node of previous elment
                EL[(i-1)*p+j-1, 0] = EL[(i-1)*p+j-2, 1]
                # third node of previous elment
                EL[(i-1)*p+j-1, 3] = EL[(i-1)*p+j-2, 2]
                EL[(i-1)*p+j-1, 1] = EL[(i-1)*p+j-1, 0] + 1
                EL[(i-1)*p+j-1, 2] = EL[(i-1)*p+j-1, 3] + 1

    """
    element_type: 2DTR3N
        2D: 2 dimensions
        TR: triangular
        3N: 3 nodes per element
    """

    if element_type == '2DTR3N':

        NoE_tri = 2*NoE        # Number of Elements
        NpE_tri = 3            # Nodes per Element

        # Generating Elements

        EL_tri = np.zeros([NoE_tri, NpE_tri])

        for i in range(1, NoE + 1):
            # first triangular element
            EL_tri[2*(i-1), 0] = EL[i-1, 0]
            EL_tri[2*(i-1), 1] = EL[i-1, 1]
            EL_tri[2*(i-1), 2] = EL[i-1, 2]

            # second triangular element
            EL_tri[2*(i-1) + 1, 0] = EL[i-1, 0]
            EL_tri[2*(i-1) + 1, 1] = EL[i-1, 2]
            EL_tri[2*(i-1) + 1, 2] = EL[i-1, 3]

        EL = EL_tri

    return NL, EL.astype('int64')
