#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_mesh.py
# Author            : David Barcene <david.barcene@utp.ac.pa>
# Date              : 10.06.2022
# Last Modified Date: 16.10.2025
# Last Modified By  : David Barcene <david.barcene@utp.ac.pa>

import numpy as np
import matplotlib.pyplot as plt
from uniform_mesh import uniform_mesh

# plt.rcParams['text.usetex'] = True  # Uncomment to render pyplot text as LaTeX
plt.rc('figure', titlesize=18)
plt.rc('axes', titlesize=15)
plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)

print('This program plots a uniform mesh of tiangular or quadrilateral finite\n\
    elements\n')

d1 = 1
d2 = 1
p = int(input("Input the number of elements in the x direction:"))
m = int(input("Input the number of elements in the y direction:"))
choose = int(input("Input the number for the element_type you want to select:\n"
                   "1. 2DQU4N\n"
                   "2. 2DTR3N\n"
                   "elemnt_type: "))

element_type = '2DTR3N'
if choose == 1:
    element_type = '2DQU4N'
elif choose == 2:
    element_type = '2DTR3N'

NL, EL = uniform_mesh(d1, d2, p, m, element_type)


NoN = np.size(NL, 0)
NoE = np.size(EL, 0)

plt.figure(figsize=(8, 8))

if element_type == '2DQU4N':
    element_counter = 1
    for j in range(NoE):
        plt.annotate(element_counter, xy=((NL[EL[j, 0]-1, 0] +
                                           NL[EL[j, 1]-1, 0] +
                                           NL[EL[j, 2]-1, 0] +
                                           NL[EL[j, 3]-1, 0])/4,
                                          (NL[EL[j, 0]-1, 1] +
                                           NL[EL[j, 1]-1, 1] +
                                           NL[EL[j, 2]-1, 1] +
                                           NL[EL[j, 3]-1, 1])/4),
                     color='blue')
        element_counter += 1

    x0 = NL[EL[:, 0]-1, 0]
    y0 = NL[EL[:, 0]-1, 1]
    x1 = NL[EL[:, 1]-1, 0]
    y1 = NL[EL[:, 1]-1, 1]
    x2 = NL[EL[:, 2]-1, 0]
    y2 = NL[EL[:, 2]-1, 1]
    x3 = NL[EL[:, 3]-1, 0]
    y3 = NL[EL[:, 3]-1, 1]

    plt.plot(np.array([x0, x1]), np.array([y0, y1]), 'red', linewidth=1)
    plt.plot(np.array([x1, x2]), np.array([y1, y2]), 'red', linewidth=1)
    plt.plot(np.array([x2, x3]), np.array([y2, y3]), 'red', linewidth=1)
    plt.plot(np.array([x3, x0]), np.array([y3, y0]), 'red', linewidth=1)


elif element_type == '2DTR3N':
    element_counter = 1
    for j in range(NoE):
        plt.annotate(element_counter, xy=((NL[EL[j, 0]-1, 0] +
                                           NL[EL[j, 1]-1, 0] +
                                           NL[EL[j, 2]-1, 0])/3,
                                          (NL[EL[j, 0]-1, 1] +
                                           NL[EL[j, 1]-1, 1] +
                                           NL[EL[j, 2]-1, 1])/3),
                     color='blue')
        element_counter += 1

    x0 = NL[EL[:, 0]-1, 0]
    y0 = NL[EL[:, 0]-1, 1]
    x1 = NL[EL[:, 1]-1, 0]
    y1 = NL[EL[:, 1]-1, 1]
    x2 = NL[EL[:, 2]-1, 0]
    y2 = NL[EL[:, 2]-1, 1]

    plt.plot(np.array([x0, x1]), np.array([y0, y1]), 'red', linewidth=1)
    plt.plot(np.array([x1, x2]), np.array([y1, y2]), 'red', linewidth=1)
    plt.plot(np.array([x2, x0]), np.array([y2, y0]), 'red', linewidth=1)

plt.show()
