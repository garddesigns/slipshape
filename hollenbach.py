#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 16:09:05 2022

@author: tobygard
"""
import numpy as np

appendages = {'bilge keels': 52}

def hollenbach(Vs, Lpp, B, T, Cb, disp):
    
    apps = {'rudder behing stern': 0.02 * Lpp * T}
    

    
    Lfore = 3.3
    Laft = 2.7
    
    Los = Lpp + Lfore + Laft
    Lc = 0
    Lwl = Lpp + Laft
    
    Ta = T
    Tf = T
    
    
    if Los <= Lpp:
        Lc = Los
    elif Los > Lpp and Los < 1.1*Lpp:
        Lc = Lpp + 2/3 * (Los - Lpp)
    else:
        Lc = 1.0667 * Lpp
        
    """Consider twin screw propellors in future"""
    
    s0 = -0.6837
    s1 = 0.2771
    s2 = 0.6542
    s3 = 0.6422
    s4 = 0.0075
    s5 = 0.0275
    s6 = -0.0045
    s7 = -0.4798
    s8 = 0.0376
    
    """Preliminary guess for propellor size"""
    
    D = 0.5962 * T
        
    k = s0 + s1 * (Los/Lwl) + s2 * (Lwl/Lpp) + s3 * Cb + s4 * (Lpp/B) + s5 * (B/T) + s6 * (Lpp/T) + s7 * ((Ta-Tf)/Lpp) + s8 * (D/T) 
    S = k * Lpp * (B + 2*T)
    
    """Frictional resistance"""
    
    Re = Vs * 0.5144 * Lc / 1.1892e-6
    
    Cf = 0.075 / (np.log10(Re)-2)**2
    
    """Residuary Resistance"""
    
    b11 = -0.57424
    b12 = 13.3893
    b13 = 90.596
    b21 = 4.6614
    b22 = -39.721
    b23 = -351.483
    b31 = -1.14215
    b32 = -12.3296
    b33 = 459.254
    
    a1 = 0.3382
    a2 = -0.8086
    a3 = -6.0258
    a4 = -3.5632
    a5 = 9.4405
    a6 = 0.0146
    a7 = 0
    a8 = 0
    a9 = 0
    a10 = 0
    
    d1 = 0.854
    d2 = -1.228
    d3 = 0.497
    e1 = 2.1701
    e2 = -0.1602
    
    Fn = Vs * 0.5144 / np.sqrt(9.81 * Lc)
    
    Cr_mean = b11 + b12 * Fn + b13 * Fn**2 + (b21 + b22 * Fn + b23 * Fn**2) * Cb + (b31 + b32 * Fn + b33 * Fn**2) * Cb**2
    
    Fr_crit = d1 + d2*Cb + d3*Cb**2
    
    c1 = (Fn/Fr_crit)
    
    kL = e1 * (Lpp)**e2
    
    if Fn < Fr_crit:
        kFr = 1
    else:
        kFr = (Fn/Fr_crit)**c1
    
    if B/T < 1.99:
        kBT = 1.99**a1
    else:
        kBT = (B/T)**a1
        
    if (Lpp/B) <= 7.11:
        kLB = (Lpp/B)**a2
    else:
        kLB = 7.11**a2
        
    if Los/Lwl <= 1.05:
        kLL = (Los/Lwl)**a3
    else:
        kLL = 1.05**a3
        
    if Lwl/Lpp <= 1.06:
        kAO = (Lwl/Lpp)**a4
    else:
        kAO = 1.05**a4

    kTr = (1 + (Ta-Tf)/Lpp)**a5
    
    if D/Ta < 0.43:
        kPr = 0.43**a6
    elif D/Ta >= 0.43 and D/Ta <= 0.84:
        kPr = (D/Ta)**a6
    else:
        kPr = 0.84**a6
        
    N_rudders = 1
    N_brackets = 2
    N_bossings = 2
    N_thrusters = 2
        
    CrBT = Cr_mean * kFr * kL * kBT * kLB * kLL * kAO * kTr * kPr * N_rudders**a7 * N_brackets**a8 * N_bossings**a9 * N_thrusters**a10
    
    Cr = CrBT * B * T / (10*S)
    
    """Correlation Allowance"""
    
    if Lpp < 175:
        Ca = (0.35-0.002*Lpp)*10**(-3)
    else:
        Ca = 0
        
    """Appendage Resistance"""
    
    k2_dict = {'rudder behind skeg': 0.35, 
               'rudder behing stern': 0.5, 
               'twin-screw rudder (slender)': 1.5, 
               'twin-screw rudder (thick)': 2.5, 
               'shaft brackets': 3.0,
               'skeg': 0.75,
               'strut bossing': 2.5,
               'hull bossings': 1.0,
               'shafts': 2.5,
               'stabiliser fins': 1.8,
               'dome': 1.7,
               'bilge keels': 0.4}
    
    sum_Sapp = 0
    k2_Sapp = 0
    for apd in apps:
        k2_Sapp += ( k2_dict[apd]) * apps[apd]
        sum_Sapp += apps[apd]
        
    k2eq = k2_Sapp/sum_Sapp
    
    Rapp = 0.5 * 1025 * (Vs*0.5144)**2 * sum_Sapp * (1 + k2eq) * Cf
    Capp = Rapp/ (0.5 * 1025 * (Vs*0.5144)**2 * S)
    
    """Environmental Resistance (Air, Wind, Wave)"""
    
    Cda = 0.8
    Avs = 383.76 #Modify!!!!!!!!
    
    Caas = Cda * 1.23 * Avs / (1025 * S )
    Cenv = Caas
    
    """Total Resistance - Only Mean fot the moment"""
    
    Ct_mean = Cf + Cr + Ca + Capp + Cenv
    
    Rt = 0.5 * 1025 * S * (Vs*0.5144)**2 * Ct_mean
    
    return Rt