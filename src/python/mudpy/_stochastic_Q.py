#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 15:13:25 2021

@author: tnye
"""

plt.figure(figsize=(6,4))
plt.loglog(f,AP)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Subfault Spectra')
# plt.ylim(8*10**-9,1.5*10**-7)
# plt.ylim(0.04,1)
plt.ylim(5.5*10**-6, 6*10**-2)
plt.show()
plt.savefig('/Users/tnye/tsuquakes/plots/path/spectra_largeQ.png', dpi=300)