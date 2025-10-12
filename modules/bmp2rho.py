#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : bmp2rho.py
# Author            : David Barcene <david.barcene@utp.ac.pa>
# Date              : 12.07.2022
# Last Modified Date: 20.07.2022
# Last Modified By  : David Barcene <david.barcene@utp.ac.pa>

import numpy as np
import cv2

img_path = './bitmap.bmp'
img = cv2.imread(img_path, 0)

img = (255.0 - img)/255.0
img = img.reshape(img.shape[0]*img.shape[1])

for i in range(len(img)):
    if img[i] != 0:
        img[i] = 1
