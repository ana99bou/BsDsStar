#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:15:26 2024

@author: anastasiaboushmelev
"""

import numpy as np
import pandas as pd
#from bokeh.plotting import figure, show, output_file
import matplotlib.pyplot as plt
import h5py
import matplotlib.pyplot as plt


f = h5py.File("BsDsStar_M2.h5", "r")
bsn0=f["/CHARM_PT_SEQ_SM10.36_s0.025/c0.340/dT26/final_state_GZ/operator_GammaZGamma5/n2_1/0_0_1/forward/data"]
bsn0b=f["/CHARM_PT_SEQ_SM10.36_s0.025/c0.340/dT26/final_state_GZ/operator_GammaZGamma5/n2_1/0_0_1/backward/data"]
ts=64
dt=26
mir = np.zeros(shape=(889, 33))
for k in range(12):
    for j in range(33):
        mir[k][j]= np.mean((np.real(bsn0[k+800, :, j]) + np.real(bsn0b[k+800, :, dt-j]))) / 2
    #mir[k][0]=(np.real(bsn0[k][0][0])+np.real(bsn0[k][1][0]))/2
   # for j in range(32):
   #     mir[k][j+1]=(np.real(bsn0[k][0][j+1])+np.real(bsn0[k][0][64-1-j])+np.real(bsn0[k][1][j+1])+np.real(bsn0[k][1][64-1-j]))/4
        
for i in range(889):
    plt.plot(range(33),mir[i],'x')
plt.yscale("log")

#########outliner: [811] bzw 812