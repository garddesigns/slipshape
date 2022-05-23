#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 17:52:39 2021

@author: tobygard
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib import cm


"""Wageningen from Ship Resitance and Propulsion"""

def KT(BAR, Z):
        KT_Wag = [[0.00880496, 0, 0, 0, 0],
               [-0.20455400, 1, 0, 0, 0],
               [0.16635100, 0, 1, 0, 0],
               [0.15811400, 0, 2, 0, 0],
               [-0.14758100, 2, 0, 1, 0], 
               [-0.48149700, 1, 1, 1, 0],
               [0.41543700, 0, 2, 1, 0],
               [0.01440430, 0, 0, 0, 1],
               [-0.05300540, 2, 0, 0, 1],
               [0.01434810, 0, 1, 0, 1],
               [0.06068260, 1, 1, 0, 1],
               [-0.01258940, 0, 0, 1, 1],
               [0.010966890, 1, 0, 1, 1],
               [-0.13369800, 0, 3, 0, 0],
               [0.00638407, 0, 6, 0, 0],
               [-0.00132718, 2, 6, 0, 0],
               [0.16849600, 3, 0, 1, 0],
               [-0.05072140, 0, 0, 2, 0],
               [0.08545590, 2, 0, 2, 0],
               [-0.05044750, 3, 0, 2, 0],
               [0.01046500, 1, 6, 2, 0],
               [-0.00648272, 2, 6, 2, 0],
               [-0.00841728, 0, 3, 0, 1],
               [0.01684240, 1, 3, 0, 1],
               [-0.00102296, 3, 3, 0, 1],
               [-0.03177910, 0, 3, 1, 1],
               [0.01860400, 1, 0, 2, 1],
               [-0.00410798, 0, 2, 2, 1],
               [-0.000606848, 0, 0, 0, 2],
               [-0.004981900, 1, 0, 0, 2],
               [0.002598300, 2, 0, 0, 2],
               [-0.000560528, 3, 0, 0, 2],
               [-0.001636520, 1, 2, 0, 2],
               [-0.000328787, 1, 6, 0, 2],
               [0.000116502, 2, 6, 0, 2],
               [0.000690904, 0, 0, 1, 2],
               [0.004217490, 0, 3, 1, 2],
               [0.0000565229, 3, 6, 1, 2],
               [-0.001465640, 0, 3, 2, 2]]
    

        PD = np.linspace(0.5 ,1.4, 10)
        J = np.linspace(0, 1.6, 30)
        
        eta_max = 0
        propellor = {}
        JKQ_funcs = []
 
        for pdr in PD:
            pdr = np.round (pdr, 1)
            KT = []
            for j in J: 
                KTj = 0
                for kt in KT_Wag:
                    KTj = KTj + kt[0] * j**kt[1] * pdr**kt[2] * BAR**kt[3] * Z**kt[4]
                KT.append(KTj)
            J1 = interp1d(KT, J)
            
            JKQ_funcs.append([pdr, J1])
        
        return JKQ_funcs
    
"""KT/KQ function defined for P/D between 0.5 and 1.4"""    
        
JKT40 = KT(0.4, 4)
JKT55 = KT(0.55, 4)
JKT70 = KT(0.7, 4)


def get_pdr(KT, BAR, J0):
    Js = []
    PDR = []
    eta_max = 0
    if(BAR == 0.4):
        for PD in JKT40:
            try:
                J = float(PD[1](KT))
                Js.append(J)
                PDR.append(PD[0])
            except:
                pass
    elif(BAR == 0.55):
        for PD in JKT55:
            try:
                J = float(PD[1](KT))
                Js.append(J)
                PDR.append(PD[0])
            except:
                pass
    elif(BAR == 0.7):
        for PD in JKT70:
            try:
                J = float(PD[1](KT))
                Js.append(J)
                PDR.append(PD[0])
            except:
                pass
        
    PD_func = interp1d(Js, PDR)
    pdr = PD_func(J0)
    
    return pdr

def KQ(BAR, J, Z, PD):
    
    KQ_Wag = [[0.00379368, 0, 0, 0, 0],
              [0.00886523, 2, 0, 0, 0],
              [-0.032241, 1, 1, 0, 0],
              [0.00344778, 0, 2, 0, 0],
              [-0.0408811, 0, 1, 1, 0],
              [-0.108009, 1, 1, 1, 0],
              [-0.0885381, 2, 1, 1, 0],
              [0.188561, 0, 2, 1, 0],
              [-0.00370871, 1, 0, 0, 1],
              [0.00513696, 0, 1, 0, 1],
              [0.0209449, 1, 1, 0, 1],
              [0.00474319, 2, 1, 0, 1],
              [-0.00723408, 2, 0, 1, 1],
              [0.00438388, 1, 1, 1, 1],
              [-0.0269403, 0, 2, 1, 1],
              [0.0558082, 3, 0, 1, 0],
              [0.0161886, 0, 3, 1, 0],
              [0.00318086, 1, 3, 1, 0],
              [0.015896, 0, 0, 2, 0],
              [0.0471729, 1, 0, 2, 0],
              [0.0196283, 3, 0, 2, 0],
              [-0.0502782, 0, 1, 2, 0],
              [-0.030055, 3, 1, 2, 0],
              [0.0417122, 2, 2, 2, 0],
              [-0.0397722, 0, 3, 2, 0],
              [-0.00350024, 0, 6, 2, 0],
              [-0.0106854, 3, 0, 0, 1],
              [0.00110903, 3, 3, 0, 1],
              [-0.000313912, 0, 6, 0, 1],
              [0.0035985, 3, 0, 1, 1],
              [-0.00142121, 0, 6, 1, 1],
              [-0.00383637, 1, 0, 2, 1],
              [0.0126803, 0, 2, 2, 1],
              [-0.00318278, 2, 3, 2, 1],
              [0.00334268, 0, 6, 2, 1],
              [-0.00183491, 1, 1, 0, 2],
              [0.000112451, 3, 2, 0, 2],
              [-0.0000297228, 3, 6, 0, 2],
              [0.000269551, 1, 0, 1, 2],
              [0.00083265, 2, 0, 1, 2],
              [0.00155334, 0, 2, 1, 2],
              [0.000302683, 0, 6, 1, 2],
              [-0.0001843, 0, 0, 2, 2],
              [-0.000425399, 0, 3, 2, 2],
              [0.0000869243, 3, 3, 2, 2],
              [-0.0004659, 0, 6, 2, 2],
              [0.0000554194, 1, 6, 2, 2]
              ]
    
    KQtot = 0
    for kq in KQ_Wag:
        KQtot = KQtot + kq[0] * J**kq[1] * PD**kq[2] * BAR**kq[3] * Z**kq[4]
        
    return KQtot
            
def prop(KT, J, Z, BAR):
    PD = get_pdr(KT, BAR, J)
    KQ0 = KQ(BAR, J, Z, PD)
    eta0 = J * KT / ( 2*np.pi*KQ0)    
    
    prop = {
        'Z': Z,
        'BAR': BAR,
        'P/D': PD,
        'J': J,
        'KT': KT,
        'KQ': KQ0,
        'eta0': eta0
        }
    
    return prop 


def propellor_design(Vs, R, L, B, T, Cb, Cp, LCB ):
    
    Z = 4
    """D limits: [0.6-0.8]*T"""
    D = T * np.linspace(0.8, 0.6, 5)
    #D = [0.6*T]
 
    """Get Estimate"""
    BAR = 0.4
    
    wt = 0.5*Cb - 0.05

    Va = Vs*(1-wt)*0.5144
    
    
    n = np.linspace(75, 200, 7)
    """Diameter for loop"""
    
    X, Y = np.meshgrid(D, n)
        
    Pd_3D = []
    
    opt_prop = {}
    Pd = 10**10
    for rpm in n:
        """Find source"""
        Pdi = []
        for d in D:
            t = (0.25014 * (B/L)**0.28956 * (np.sqrt(B*T)/d)**0.2624)/((1 - Cp + 0.0225*LCB)**0.01762)
            thrust = R/(1-t)
            KT = thrust/(1025 * (rpm/60)**2 * d**4)
            J = Va/((rpm/60)*d)
            try:
                nr = rpm/60
                propellor = prop(KT, J, Z, BAR)
                
                
                eta_0 =propellor['eta0']
                eta_h = (1-t)/(1-wt)
                eta_r = 0.9922 - 0.05908 * BAR + 0.07424 * (Cp - 0.0225*LCB)
                Pd1 = 2 * np.pi * (rpm/60) * propellor['KQ'] * 1025 * (rpm/60)**2 * d**5
                """Pd1 = R*Vs*0.5144/(eta_0 * eta_r * eta_h)"""
                
                if(Pd1<Pd):
                    Pd = Pd1
                    opt_prop = {
                        'rpm': rpm,
                        'D': d,
                        'eta0': eta_0,
                        'etar': eta_r,
                        'etah': eta_h,
                        'Va': Va,
                        'KT':KT,
                        'P/D': propellor['P/D'],
                        'Pd': Pd,
                        'L': L,
                        'B': B,
                        'T': T
                        }
            except:
                Pd1 = 'NaN'
            
            cav_check = False
            Pdi.append(Pd1)
        Pd_3D.append(Pdi)
    
    """Pd_3D = np.array(Pd_3D)
    
    fig = plt.figure()
    ax = plt.axes()
    cs = ax.contourf(X,Y, Pd_3D, 20, cmap="RdBu")
    ax.set_xlabel('Diameter(m)')
    ax.set_ylabel('Revolutions/minutes')
    cbar = fig.colorbar(cs)
    plt.savefig(str(np.round(R))+'.png')
    plt.show()"""

    return Pd


def opt_rpm(Vs, Pe, wt, t, BAR, Z, D, T):
    vs = 0.5144 * Vs
    
    Thrust = (Pe/vs)/(1-t)
    va = vs * (1-wt)
    
    h = T - 0.5*D
    
    def KT(n):
        return Thrust/(1025* n**2 * D**4)
    
    def J(n):
        return va/(n * D)

    rpm = np.linspace(80, 200, 13)
    #rpm = np.array([102, 114, 126])
    n = rpm/60
    
    eta0 = []
    min_PD = 10e9
    opt_prop = {}
    
    for i in n:
        try:
            kt = KT(i)
            j = J(i)
            
            propellor = prop(kt, j, Z, BAR)
            pdr = propellor['P/D']
            eta0.append(propellor['eta0'])
            
            req_bar = cavitation(Thrust, va, i, D, pdr, h)
            
            
            if 0.4< req_bar < 0.55:
                req_bar1 = req_bar
                cav_check = True
            elif 0.55<= req_bar< 0.7:
                propellor = prop(kt, j, Z, 0.55)
                req_bar1 = cavitation(Thrust, va, i, D, propellor['P/D'], h)
                cav_check = True
            elif req_bar > 0.7:
                pp = prop(kt, j, Z, 0.7)
                req_bar1 = cavitation(Thrust, va, i, D, propellor['P/D'], h)
                cav_check = True
            else:
                cav_check = False
            
            #print(f'n = {i:2g}rps, J = {j:3g}, KT = {kt:3g}, P/D = {pdr:3g}, eta0 = {propellor["eta0"]:3g}')
            
            pd = 2 * np.pi * i * propellor['KQ'] * 1025 * i**2 * D**5
             
            if pd < min_PD:
                min_PD = pd
                opt_prop = {
                    'Va': va,
                    'n': i,
                    'D': D,
                    'KT': kt,
                    'J': j,
                    'P/D': float(pdr),
                    'Thrust': Thrust,
                    'BAR': BAR,
                    'Required BAR':req_bar1,
                    'h': h
                    }
        except:
            idx = np.where(n == i)
            n = np.delete(n, idx)
    
    """fig, ax = plt.subplots()
    ax.plot(n, eta0)
    fig.show()"""
    
    return min_PD, opt_prop #return max eta and rpm

def cavitation(thrust, Va, n, D, PD, h):
     
    Vr = np.sqrt(Va**2+(0.7*np.pi*n*D)**2)
    
    sigma = (1025*9.81*h + 101*10**3 - 3*10**3)/(0.5*1025*Vr**2)
    
    # Merchant Ships
    tau_c = 0.28*(sigma - 0.03)**0.57
    
    Ap = thrust/(0.5*1025*Vr**2*tau_c)
    Ad = Ap/(1.067-0.229*PD)
    BAR = Ad/(np.pi*D**2/4)
    return BAR

def prel_screen(Vs, Pe, wt, BAR, Z, L, B, T, Cp, LCB):
    D = np.linspace(0.6*T, 0.8*T, 5)
    
    #fig, ax = plt.subplots()
    
    min_PD = 10e9
    opt_prop = {}
    
    for d in D:
        t = (0.25014 * (B/L)**0.28956 * (np.sqrt(B*T)/d)**0.2624)/((1 - Cp + 0.0225*LCB)**0.01762)
        PD, pp = opt_rpm(Vs, Pe, wt, t, BAR, Z, d, T)
            
        if PD < min_PD:
            min_PD = PD
            opt_prop = pp
            
        #ax.plot(n, eta0, label=d)
        
    """fig.legend(title="Propellor Diametors")
    ax.set_ylabel(r'$\eta_0$')
    ax.set_xlabel('n(rps)')
    ax.grid()
    fig.show()"""
    
    return min_PD, opt_prop
    


def propellor_design(Vs, Pe, L, B, T, Cb, Cp, LCB):
    
    wt = 0.5*Cb - 0.05
    BAR = 0.4
    Z= 4
    
    PD, opt_prop = prel_screen(Vs, Pe, wt, BAR, Z, L, B, T, Cp, LCB)
    
    
    return PD
