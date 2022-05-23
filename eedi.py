#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 16:29:36 2021

@author: tobygard
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def EEDI(Pi, C, V, fuel):
    
    Cf = {
        "diesel":3.206,
        "LFO": 3.151,
        "HFO": 3.114,
        "LPG": 3.000,
        "LNG": 3.030,
        "methanol": 1.375,
        "ethanol": 1.913
        }
    Pme = (0.75 * Pi)/10**3
    
    if Pme > 10000:
        Pae = 0.025 * Pme + 250
    else:
        Pae = 0.05 * Pme
    
    """Get Data"""
    
    SFC_ME = 171
    SFC_AE = 215
    
    eedi = (Pme * Cf[fuel] * SFC_ME + Pae * Cf[fuel] * SFC_AE)/ (C * V)
    
    return eedi

eedi_refs = {
    'bulker': {'a':961.79, 'c': 0.477},
    'tanker': {'a':1218.80, 'c': 0.488},
    'cargo': {'a':107.48, 'c': 0.216}
    }

def ref_eedi(C, ship_type):
    return eedi_refs[ship_type]['a']* C **(-eedi_refs[ship_type]['c'])

def EEDI_plot(min_eedi, max_eedi, DWT, ship_type, L_min, L_max, B_min, B_max, T_min, T_max, Cb_min, Cb_max):
    
    Min = "(" + str(int(L_min)) + ", " + str(int(B_min)) + ", " + str(int(T_min)) + ", " + str(round(Cb_min, 2)) + ")"
    Max = "(" + str(int(L_max)) + ", " + str(int(B_max)) + ", " + str(int(T_max)) + ", " + str(round(Cb_max, 2)) + ")"
    
    fig, ax = plt.subplots()
    plt.ylim(ymin=0, ymax=17)
    
    capacities = np.linspace(100, 1.4 * DWT, 40)
    ax.plot(capacities, ref_eedi(capacities, ship_type))
    # hullsdf['EEDI'].max()
    # 
    ax.plot(DWT, float(max_eedi), 'b.')
    ax.plot(DWT, float(min_eedi), 'r.')
    ax.annotate(Max, (DWT, float(max_eedi)))
    ax.annotate(Min, (DWT, float(min_eedi)))
    fig.savefig('static/images/eedi.png')
    
def time_eedi(DWT, year):
    Tier1 = 0.1 # 2015  - 2019
    Tier2 = 0.15 # 2019 - 2025
    Tier3 = 0.3 # 2025 - onwards
    
    t0 = ref_eedi(DWT, 'bulker')
    t1 = (1-Tier1) * t0
    t2 = (1-Tier2) * t0
    t3 = (1-Tier3) * t0
    
    years = [2013, 2015, 2019, 2025]
    eedis = [t0, t1, t2, t3]
    
    eedi = interp1d(years, eedis, fill_value="extrapolate")
    
    
    return eedi(year)
    
def opt_hull_eedi(eedi, hulls):
    below_eedi_hulls = hulls[hulls['EEDI']<eedi]
    
    max_idx = below_eedi_hulls['EEDI'].idxmax()
    L_max, B_max, T_max, Cb_max = below_eedi_hulls['L'].mean(), below_eedi_hulls['B'].mean(), below_eedi_hulls['T'].mean(), below_eedi_hulls['Cb'].mean()
    
    
    return L_max, B_max, T_max, Cb_max

def percentage_eedi(eedi, years, hulls, DWT):
     
    L = []
    B = []
    T = []
    Cb =[]
    
    for i in eedi:
        new_db = hulls[hulls['EEDI']<i]
        L.append(new_db['L'].to_numpy())
        B.append(new_db['B'].to_numpy())
        T.append(new_db['T'].to_numpy())
        Cb.append(new_db['Cb'].to_numpy())
   
    dlt_idx = [] 
    for j in range(len(eedi)):
        if len(L[j]) == 0:
            dlt_idx.append(j)
        
    L = np.delete(L, dlt_idx)
    eedi = np.delete(eedi, dlt_idx)
    years = np.delete(years, dlt_idx)
    Cb = np.delete(Cb, dlt_idx)
    T = np.delete(T, dlt_idx)
    B = np.delete(B, dlt_idx)
    
    fig, ax = plt.subplots()

    ax.boxplot(L, labels=years)
    ax.set_xlabel('Year')
    ax.set_ylabel('Length(m)')
    
    fig.savefig(f'static/images/EEDI/L DWT={DWT}.png')
    
    fig, ax = plt.subplots()

    ax.boxplot(B, labels=years)
    ax.set_xlabel('Year')
    ax.set_ylabel('Breadth(m)')
    
    fig.savefig(f"static/images/EEDI/B DWT={DWT}.png")
    
    fig, ax = plt.subplots()

    ax.boxplot(T, labels=years)
    ax.set_xlabel('Year')
    ax.set_ylabel('Draught(m)')
    
    fig.savefig(f"static/images/EEDI/T DWT={DWT}.png")
    
    fig, ax = plt.subplots()

    ax.boxplot(Cb, labels=years)
    ax.set_xlabel('Year')
    ax.set_ylabel('Block Coefficient')
    
    fig.savefig(f"static/images/EEDI/Cb DWT={DWT}.png")

        
        
        