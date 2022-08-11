#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 14:50:02 2022

@author: dmelgarm
"""
import numpy as np
from matplotlib import pyplot as plt

f = np.genfromtxt('/Users/dmelgarm/FakeQuakes/tohoku_simple/output/ruptures/minson.000000.rupt')

plt.figure()
plt.subplot(131)
plt.scatter(f[:,1],f[:,2],c=f[:,9],cmap='jet')
plt.colorbar(label='Slip (m)')

plt.subplot(132)
plt.scatter(f[:,1],f[:,2],c=f[:,7],cmap='magma')
plt.colorbar(label='Rise time (s)')

plt.subplot(133)
plt.scatter(f[:,1],f[:,2],c=f[:,12],cmap='hot')
plt.colorbar(label='Onset time (s)')
