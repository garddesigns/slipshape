# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 14:16:58 2022

@author: tg4g19
"""

import matplotlib.pyplot as plt
import numpy as np

"""----------------------- DISPLACEMENT METHODS ----------------------------"""

dictionnary = {
    "cargo": {
        "limits": {"DWT": [[5000, 15000]]},
        "DWT/Disp": [[0.65, 0.8]],
        "Wst/Wl": [[0.55, 0.64]],
        "Wot/Wl": [[0.19, 0.33]],
        "Wm/Wl": [[0.11, 0.22]]
        },
    "coaster": {
        "limits": {"GRT": [[499, 999]]},
        "DWT/Disp": [[0.70, 0.75]],
        "Wst/Wl": [[0.57, 0.62]],
        "Wot/Wl": [[0.30, 0.33]],
        "Wm/Wl": [[0.09, 0.12]]
        },
    "bulker": {
        "limits": {"DWT": [[20000, 50000], [50000, 200000]]},
        "DWT/Disp": [[0.74, 0.85],[0.8,0.87]],
        "Wst/Wl": [[0.68, 0.79],[0.78, 0.85]],
        "Wot/Wl": [[0.10, 0.17],[0.06, 0.13]],
        "Wm/Wl": [[0.12, 0.16], [0.08, 0.14]]
        },
    "tanker": {
        "limits": {"DWT": [[25000, 120000], [200000, 500000]]},
        "DWT/Disp": [[0.78, 0.86],[0.83, 0.88]],
        "Wst/Wl": [[0.73, 0.83],[0.75, 0.88]],
        "Wot/Wl": [[0.05, 0.12],[0.09, 0.13]],
        "Wm/Wl": [[0.11, 0.16], [0.09, 0.16]]
        },
    "containership": {
        "limits": {"DWT": [[10000, 15000], [15000, 165000]]},
        "DWT/Disp": [[0.65, 0.74],[0.65, 0.76]],
        "Wst/Wl": [[0.58, 0.71],[0.62, 0.72]],
        "Wot/Wl": [[0.15, 0.20],[0.14, 0.20]],
        "Wm/Wl": [[0.09, 0.22], [0.15, 0.18]]
        },
    "ro-ro": {
        "limits": {"DWT": [[0, 16000]], "L":[[80, 500]]},
        "DWT/Disp": [[0.50, 0.60]],
        "Wst/Wl": [[0.68, 0.78]],
        "Wot/Wl": [[0.12, 0.19]],
        "Wm/Wl": [[0.10, 0.20]]
        },
    "reefer": {
        "limits": {"Displacement": [[8495, 14158]]},
        "DWT/Disp": [[0.45, 0.55]],
        "Wst/Wl": [[0.51, 0.62]],
        "Wot/Wl": [[0.21, 0.28]],
        "Wm/Wl": [[0.15, 0.26]]
        },
    "cruise": {
        "limits": {"L": [[200, 360]]},
        "DWT/Disp": [[0.23, 0.34]],
        "Wst/Wl": [[0.52, 0.56]],
        "Wot/Wl": [[0.30, 0.34]],
        "Wm/Wl": [[0.15, 0.20]]
        },
    "trawler": {
        "limits": {"L": [[44, 82]]},
        "DWT/Disp": [[0.30, 0.58]],
        "Wst/Wl": [[0.42, 0.46]],
        "Wot/Wl": [[0.36, 0.40]],
        "Wm/Wl": [[0.15, 0.20]]
        },
    "tug": {
        "limits": {"Pb": [500000, 3000000]},
        "DWT/Disp": [0.20, 0.40],
        "Wst/Wl": [0.42, 0.56],
        "Wot/Wl": [0.17, 0.21],
        "Wm/Wl": [0.38, 0.43]
        }
    
    }

def displacement(DWT, n, ship_type="cargo"):
    """n = number of displacement ratio iterations"""
    displacements = []
    if(list(dictionnary[ship_type]["limits"].keys()).count('DWT') > 0):
        for i in range(len(dictionnary[ship_type]["limits"]["DWT"])):
            if(DWT > dictionnary[ship_type]["limits"]["DWT"][i][0] and DWT < dictionnary[ship_type]["limits"]["DWT"][i][1]):
               disp_ratios = np.linspace(dictionnary[ship_type]["DWT/Disp"][i][0], dictionnary[ship_type]["DWT/Disp"][i][1], n)
               displacements = DWT / disp_ratios
    else:
        disp_ratios = np.linspace(dictionnary[ship_type]["DWT/Disp"][0][0], dictionnary[ship_type]["DWT/Disp"][0][1], n)
        displacements = DWT / disp_ratios
        
    if(list(dictionnary[ship_type]["limits"].keys()).count('Displacement') > 0):
        filtered = filter(lambda disp: disp > dictionnary[ship_type]["limits"]["Displacement"][0][0] and disp < dictionnary[ship_type]["limits"]["Displacement"][0][1], displacements)
        displacements = filtered
        
    displacements = list(displacements)
    
    if displacements:
        return displacements
    else:
        return "DWT out of range for selected vessel type" 
    
"""-------------------------- ITERATIVE METHODS ----------------------------"""

dictionary2 = {
    "cargo":{
        "limits": {'Fn': [[0.25, 0.5], [0, 0.25]]},
        "Cb": [[0.56, 0.64], [0.65, 0.73]],
        "Lpp/B": [[5.7, 7.8], [4.8, 8.5]],
        "B/T": [[2.2, 2.6], [2.1, 2.3]],
        "Lpp/Disp": [[5.6, 5.9], [5.2, 5.4]]
        },
    "coaster":{
        "Cb": [[0.58, 0.72]],
        "Lpp/B": [[4.5, 5.5]],
        "B/T": [[2.5, 2.7]],
        "Lpp/Disp": [[4.2, 4.8]]
        },
    "tug":{
        "Cb": [[0.5, 0.58]],
        "Lpp/B": [[3.8, 4.5]],
        "B/T": [[2.4, 2.6]],
        "Lpp/Disp": [[4.0, 4.6]]
        },
    "tanker":{
        "limits": {'Fn': [[0, 0.2], [0.14, 0.25]]},
        "Cb": [[0.82, 0.88], [0.78, 0.86]],
        "Lpp/B": [[5.1, 6.8], [5.0, 6.5]],
        "B/T": [[2.4, 3.2], [2.2, 2.9]],
        "Lpp/Disp": [[4.5, 5.6], [6.1, 6.5]]
        },
    "reefer":{
        "Cb": [[0.57, 0.59]],
        "Lpp/B": [[6.7, 7.2]],
        "B/T": [[2.8, 3.0]],
        "Lpp/Disp": [[6.1, 6.5]]
        },
    "bulker":{
        "Cb": [[0.65, 0.92]],
        "Lpp/B": [[5.0, 7.1]],
        "B/T": [[2.1, 3.2]],
        "Lpp/Disp": [[4.7, 5.6]]
        }
    }


def preliminary_dimensions(disp, V, n, ship_type='cargo'):
    
    
    """Length"""
    ships = []
    for d in disp:
        for i in range(len(dictionary2[ship_type]["Cb"])):
            L = np.linspace(dictionary2[ship_type]["Lpp/Disp"][i][0], dictionary2[ship_type]["Lpp/Disp"][i][1], n)*d**(1/3)
            Fn = V*0.5144/np.sqrt(9.81*L)
            if("limits" in dictionary2[ship_type]):
                 """Fn = np.linspace(dictionary2[ship_type]["limits"]["Fn"][i][0], dictionary2[ship_type]["limits"]["Fn"][i][1], n)
                 L = (V*0.5144)**2/(9.81*Fn**2)
                 for j in L:
                     print(j/d**(1/3), dictionary2[ship_type]["Lpp/Disp"][i][0], dictionary2[ship_type]["Lpp/Disp"][i][1])
                     if j/d**(1/3) < dictionary2[ship_type]["Lpp/Disp"][i][0] or j/d**(1/3) > dictionary2[ship_type]["Lpp/Disp"][i][1]:
                         L = np.delete(L, np.where(L == j )[0][0])"""
                 for j in Fn:
                    if j < dictionary2[ship_type]["limits"]["Fn"][i][0] or j > dictionary2[ship_type]["limits"]["Fn"][i][1]:
                        Fn = np.delete(Fn, np.where(Fn == j)[0][0])
                        
                        
            L = (V*0.5144)**2/(9.81*Fn**2)
            
            
            if L.size == 0:
                continue
            
            Cb = np.linspace(dictionary2[ship_type]["Cb"][i][0], dictionary2[ship_type]["Cb"][i][1], n)
            
            
            for l in L:
                B = l/np.linspace(dictionary2[ship_type]["Lpp/B"][i][0], dictionary2[ship_type]["Lpp/B"][i][1], n)
                for b in B:
                    for cb in Cb:
                        T = (d/1.025)/(l*b*cb)
                        if b/T < dictionary2[ship_type]["B/T"][i][0]  or b/T > dictionary2[ship_type]["B/T"][i][1]:
                            continue
                        else:
                            ships.append([d, l, b, T, cb])
                        #ships.append([d, l, b, T, cb])
                                         
    return ships


"""----------------------------- LENGTH METHODS ----------------------------"""

def schneekluth_length(disp, V, Cb):
    
    """Conditions of use:
        -Displacemeent > 1000 tonnes
        -0.16 <= Fn <= 0.32
        -0.48 <= Cb <= 0.85
        """
    
    Fn = np.linspace(0.16, 0.32, 200)
    L0 = (V*0.5144)**2/(9.81*Fn**2)
    
    def L(Fn):
        arr = []
        for fn in Fn:
            if((0.145/fn)==Cb):
                arr.append(disp**0.3 * V**0.3 * 3.2)
            else:
                C = 3.2*(Cb + 0.5)/((0.145/fn) + 0.5)
                arr.append(disp**0.3 * V**0.3 * C)
                
        return arr
        
    L1 = L(Fn) 
    index = np.argwhere(np.diff(np.sign(L1-L0))).flatten()    
  
    
    fig, ax = plt.subplots()
    ax.plot(Fn, L0, '-')
    ax.plot(Fn, L1, '-')
    ax.plot(Fn[index], L0[index], 'ko')
    plt.show()
    
    return Fn[index], L0[index]

def ayre_length(disp, V):
    
    Fn = np.linspace(0.1, 0.5, 200)
    Lpp = (V*0.5144)**2/(9.81*Fn**2)
    f = Lpp/((disp/1.025)**(1/3))
    g = 3.33 + 1.67 * (V/np.sqrt(Lpp))
    
    index = np.argwhere(np.diff(np.sign(f-g))).flatten() 
    
    ig, ax = plt.subplots()
    ax.plot(Lpp, f, '-')
    ax.plot(Lpp, g, '-')
    ax.plot(Lpp[index], f[index], 'ko')
    plt.show()
    
    return Lpp[index]

def posdunine_length(disp, V):
    """Unrealsistic answers DO NOT USE AT THE MOMENT"""
    return ((disp/1.025)**(1/3)*7.62*V)/((V+2)**2)


def volker_length(disp, V, ship="cargo"):
    C1 = 0
    if(ship=="cargo" or ship=="container" or ship=="bulker"):
        C1 = 3.5
    elif(ship=="reefer"):
        C1 = 3.0
    elif(ship=="fishing"):
        C1 = 2.0
        
    return (C1 + (4.5 * V * 0.5144)/np.sqrt(9.81*(disp/1.025)**(1/3)))*(disp/1.025)**(1/3)

"""---------------------------- BREADTH METHODS ----------------------------"""

def empirical_breadth(L, ship_type="cargo"):
    dictionary_b = {
        "cargo": {
            "length": [50, 200],
            "L/B": 4 + 0.015 * (L + 17)
            },
        "reefer": {
            "length": [60, 180],
            "L/B": 4 + 0.014 * (L + 11)
            },
        "container": {},
        "bulker": {
            "length": [120, 500],
            "L/B": 6
            },
        "tanker":  {
            "length": [0, 500],
            "L/B": 5.5
            },
        "LPG": {
            "length": [100, 500],
            "L/B": 5.7 + 0.002 * (L - 100)
            },
        "LNG": {
            "length": [100, 500],
            "L/B": 5.7 + 0.002 * (L - 100)
            },
        "RoRo": {
            "length": [80, 500],
            "L/B": 5.5 + 0.0036 * (L - 41)
            },
        "ROPAX": {
            "length": [80, 500],
            "L/B": 5.7 + 0.0033 * (L - 141)
            }
        }
    
    if(L < dictionary_b[ship_type]["length"][0]):
        return "Length too short for selected ship type"
    elif(L > dictionary_b[ship_type]["length"][1]):
        return "Length too large for selected ship type"
    else:
        return L/dictionary_b[ship_type]["L/B"]
    
"""------------------ BLOCK COEFFICIENT METHODS ----------------------------"""

"""Section 2.10.6 Table 2.15/2.16"""

def empirical1_Cb(Fn):
    """For single screw ships at service speed, assuming service speed and 
    trial speed are different"""
    K1 = 1.06
    K2 = 1.68
    K3 = 0
    
    return K1 - K2 * Fn - K3**2 * Fn

def empirical2_Cb(V, L):
    """Assuming "Anglo-Saxon" units = imperial units
    Only for Cargo Ships"""
    
    L = 3.28084 * L
    F = V / np.sqrt(L)
    
    K4 = 0
    K5 = 0
    
    if(F > 0.65 and F < 0.8):
        "Cargo Ships"
        K4 = 1.06
        K5 = 0.500
    elif(F > 0.89):
        "Fast Cargo Ships"
        K4 = 1.03
        K5 = 0.500
    elif(F < 0.65):
        K4 = 1.12
        K5 = 0.500
    else:
        return  "Out of applicable range"
    
    return K4 - K5 * F
        
    

"""Section 2.10.6 Part B (method 2 and 3)"""

def Schneekluth1_Cb(V,L,B):
    Fn = (V*0.5144)/(9.81*L)**(1/2)
    Cb = (0.14 * (L/B) + 20)/(Fn * 26)
    if Cb < 0.48 or Cb > 0.85:
        return "Out of Cb range"
    if Fn < 0.14 or Fn > 0.32:
        return "Out of appropriate Fn range "
    return Cb

"""------------------------- MISC ------------------------------------------"""

"""Only valid for tankers and bulkers"""
def A_v(Loa):
    FB = 10.65 + 0.0515 * Loa #Freeboard
    H = 9.11 + 0.0260 * Loa #Bridge to Deck Height
    return FB, H
    