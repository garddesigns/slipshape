 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 10:50:41 2021

@author: tobygard
"""

import pandas as pd
import numpy as np
from preliminary_dimensions import A_v

example_data = [[0.23655, 18872, 0.6492, 147.7, 24, 8.2, 1.16e9,  1.3067, 0.7675, 4400]]
data1 = pd.DataFrame(example_data, columns = ['Fn', 'Displacement', 'Cb', 'L', 'B', 'T', 'Rn', "LCB0", 'Cwp0', 'S0'])



def holtrop_resistance(disp, V, Cb, L, B, T, ab=0):
    
    apps = {'rudder behing stern': 0.02 * L * T}
    
    
    Fn = V*0.5144/np.sqrt(9.81*L)
    
    
    Rn = V*0.5144*L/(1.16*10**-6)
    
    
    """Assumptions: as A_T = 0 (false), Ta = T """
    rho_water = 1025
    rho_air = 1.225
    A_T = 0
    A_BT = 0
    h_B = 0
    Ta = 0
    """Find estiamtion method"""
    FB, H = A_v(L)
    A_front = (FB+H) * B
    """Normal Section"""
    C_stern = 0
    
    """Holtrop Resistance Paper estimations"""
    
    LCB = -(0.44*Fn-0.094)
    
    """ Method - Laboratory HSVA"""
    Cm = 1/(1+(1-Cb)**3.5)
    #print(Cm)
    Cp = (Cb/Cm)
    
    
    """For tanker, bulk carrier, genral cargo with Cp 0.56<Cp<0.87"""
    Cwp = 0.763*(Cp+0.34)
    
    
    """Wetted surface area, no bulbous bow (A_BT = 0)"""
    
    S = L * (2*T+B) * Cm**0.5 * (0.453 + 0.4425 * Cb - 0.2862 * Cm - 0.003467 * (B/T) + 0.3696 * Cwp) + 2.38 * A_T + (A_BT/Cb)
    """Waterline entrance angle"""
    
    Lr = L * ((1-Cp)+ ((0.06* Cp * LCB)/(4 * Cp- 1)))
    a = -((L/B)**0.80856 * (1-Cwp)**0.30484 * (1-Cp-0.0225*LCB)**0.6367 * (Lr/B)**0.34574 * ((100*disp/(L)**3)**0.16302))
    iE = 1 + 89* np.exp(a)
    
    
    
    """LOOK INTO THIS!!!!!!"""
    
    Tf = T
    
    
    """-----------------------------------------------RESISTANCE ESTIMATION ------------------------------------------------------------------"""
    
    """Hull form factor"""
    
    c14 = 1 + 0.011*C_stern
    k = -0.07 + 0.487118 * c14 * ((B/L)**1.06806 * (T/L)**0.46106 * (L/Lr)**0.121563 * (L**3/disp)**0.36486 * (1- Cp)**(-0.604247)) 
    
    #print(f'Form Factor : {k:.3g}')
    
    """Frictional Resistance with ITTC correlation line"""
    
    
    Cf = 0.075/(((np.log10(Rn)-2)**2))
    
    """ Assume salt water (density =1.025)"""
    
    def air_lube(Cf_lube, airlubearea):
        S_bottom = L * B * Cm**0.5 * (0.453 + 0.4425 * Cb - 0.2862 * Cm - 0.003467 * (B/T) + 0.3696 * Cwp)
        S_side = L * 2*T * Cm**0.5 * (0.453 + 0.4425 * Cb - 0.2862 * Cm - 0.003467 * (B/T) + 0.3696 * Cwp) 
        S_transom = 2.38 * A_T 
        
        S_lube = airlubearea * S_bottom
        S_nonlube = (1-airlubearea) * S_bottom + S_side + S_transom
        
        R_lube = 0.5 * rho_water * (V*0.5144)**2 * S_lube * Cf * (1-Cf_lube)
        R_nonlube = 0.5 * rho_water * (V*0.5144)**2 * S_nonlube * Cf 
        
        Rf = R_lube + R_nonlube
        
        return Rf
    
    if ab:
        Rf = air_lube(ab['Cf_lube'], ab['airlubearea'])
    else:
        Rf = 0.5 * rho_water * (V*0.5144)**2 * S * Cf
        #print(f'Frictional Resistance: {round(Rf)} N')
    
    
    """Appendage Resistance : ASSUMPTIONS TO BE MADE """
    
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
        k2_Sapp += (1 + k2_dict[apd]) * apps[apd]
        sum_Sapp += apps[apd]
        
    k2eq = k2_Sapp/sum_Sapp
    
    Rapp = 0.5 * rho_water * (V*0.5144)**2 * sum_Sapp * k2eq * Cf
    
    #print(f'Appendage Resistance: {round(Rapp)} N')

    """Wave Resistance"""
    
    if((B/L) <=0.11):
        c7 = 0.229577 * (B/L)**(1/3) 
    elif(0.11< (B/L) <= 0.25):
        c7 = (B/L) 
    elif((B/L) > 0.25):
        c7 = 0.5 - 0.0625 * (B/L) 
        
    c1 = 2223105 * (c7**3.78613) * ((T/B)**1.07961) * ((90 - iE)**(-1.37565))
    
    """No bulbous bow, Tf = T/2.572 (needs looking into)"""
    
    c3 = 0.56 * (A_BT**1.5)/(B*T*(0.31*(A_BT**0.5) + Tf - h_B))
    c2 = np.exp(-1.89*np.sqrt(c3))
        
    if((L**3/disp) <= 512):
        c15 = 1.69385*(-1) 
    elif(512 < (L**3/disp) <= 1726.91):
        c15 = 1.69385*(-1) - ((L/disp**(1/3))-8)/2.36
    elif((L**3/disp) > 1726.91):
        c15 = 0
    
    if(Cp <= 0.8):
        c16 = 8.07981 * Cp - 13.8673 * Cp**2 + 6.984388 * Cp**3
    elif(Cp > 0.8):
        c16 = 1.73014 - 0.7067*Cp
        
    if((L/B) <= 12):
        Lambda = 1.446 * Cp - 0.03 * (L/B)
    elif((L/B) > 12):
        Lambda = 1.446 * Cp - 0.36
    
    c17 = 6919.3 * Cm**(-1.3346) * (disp/L**3)**2.00977 * ((L/B)-2)**1.40692
    
    m1 = 0.0140407*(L/T) - 1.75254*(disp**(1/3)/L) - 4.79323*(B/L) - c16
    m4 = 0.4 * c15 * np.exp(-0.034 * Fn**(-3.29))
    m3 = -7.2035 * (B/L)**0.326869 * (T/B)**0.605375
    d = -0.9
    c5 = 1 - 0.8* A_T/(B * T * Cm)
    
    if(Fn <= 0.4):
        Rw = c1 * c2 * c5 * rho_water * 9.81 * disp * np.exp(m1 * ((Fn)**d) + m4 * np.cos(Lambda*(Fn**(-2))))
    elif(0.4< Fn <= 0.55):
        m4a = 0.4 * c15 * np.exp(-0.034 * 0.4**(-3.29))
        m4b = 0.4 * c15 * np.exp(-0.034 * 0.55**(-3.29))
        Rwa = c1 * c2 * c5 * rho_water * 9.81 * disp * np.exp(m1* (0.4**d) + m4a * np.cos(Lambda*(0.4**(-2))))
        Rwb = c17 * c2 * c5 * rho_water * 9.81 * disp * np.exp(m3* (0.55**d) + m4b * np.cos(Lambda*(0.55**(-2))))
        Rw = Rwa + ((20 * Fn - 8)/3) * (Rwb-Rwa)
    elif(Fn > 0.55):
        Rw = c17 * c2 * c5 * rho_water * 9.81 * disp * np.exp(m3* ((Fn)**d) + m4 * np.cos(Lambda*(Fn**(-2))))
    
    #print(f'Wave Resistance: {round(Rw)} N')
    
    """Correlation allowance resistance"""
    
    if((T/L) <= 0.04):
        c4 = (T/L)
    elif((T/L) > 0.04):
        c4 = 0.04
    
    """Correlation allowance coefficient - No Roughness Allowance as new vessel"""
    
    C_A = 0.006 * (L + 100)**(-0.16) - 0.00205 + 0.003 * np.sqrt(L/7.5) * Cb**4 * c2 * (0.04 - c4)
    
    RA = 0.5 * rho_water * (V*0.5144)**2 * C_A * (S + sum_Sapp)
    
    #print(f'Correlation Allowance: {round(RA)} N')
    
    
    
    """Air resistance"""
    
    C_DA = 0.8 
    RAA = 0.5 * rho_air * (V*0.5144)**2 * C_DA * A_front
    
    #print(f'Air Resistance: {round(RAA)} N')
    
    """Total Resistance, remove 1.1 factor"""
    Rt = Rf*(1+k) + Rw + RA + Rapp + RAA
    
    Cr = (Rt-Rf)/ (0.5 * rho_water * (V*0.5144)**2*S)
    
    #print("-------------------------------------------------------------------------------")
   
 
    return (Rt, Cp, LCB, Cwp, Cr)
