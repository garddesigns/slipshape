#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 14:00:43 2021

@author: tobygard

"""

import pandas as pd
import numpy as np 
from holtrop import holtrop_resistance
from eedi import EEDI
from eedi import EEDI_plot
from eedi import opt_hull_eedi
from eedi import time_eedi
from eedi import percentage_eedi
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import preliminary_dimensions as preldim
from propdesign import propellor_design 

def shipdesign(V, DWT, vessel_type, fuel, air_lub={}, dimensions=[], plots=[], n=7):
    
    displacement = preldim.displacement(DWT, n, vessel_type)
        
    hulls = preldim.preliminary_dimensions(displacement, V, n, vessel_type)

    hullsdf = pd.DataFrame(hulls, columns = ['Displacement', 'L', 'B', 'T', 'Cb'])
    
    hullsdf['T'] = hullsdf['Displacement']/(hullsdf['L']*hullsdf['B']*hullsdf['Cb']*1.025)
    hullsdf['D^(1/3)'] =  hullsdf['Displacement']**(1/3)
    
    """Assumed rudder is 2% of LT"""
    
    hullsdf['A_rudder'] = 0.02 * hullsdf['L'] * hullsdf['T']
    
    def Rt(df):
        return holtrop_resistance(df['Displacement'], V, df['Cb'], df['L'], df['B'], df['T'], air_lub)[0]
        
    def Cp(df):
        return holtrop_resistance(df['Displacement'], V, df['Cb'], df['L'], df['B'], df['T'], air_lub)[1]
        
    def LCB(df):
        return holtrop_resistance(df['Displacement'], V, df['Cb'], df['L'], df['B'], df['T'], air_lub)[2]
    
    def S(df):
        return holtrop_resistance(df['Displacement'], V, df['Cb'], df['L'], df['B'], df['T'], air_lub)[4]
    
    "Do ...[0], ...[1] for Cp, LCB "
    hullsdf['Rt'] = hullsdf.apply(Rt, axis = 1)
    hullsdf['LCB'] = hullsdf.apply(LCB, axis = 1)
    hullsdf['Cp'] = hullsdf.apply(Cp, axis = 1)
    hullsdf['Pe'] = hullsdf['Rt']*0.5144*V
    hullsdf['S'] = hullsdf.apply(S, axis = 1)

    
    def Pd(df):
        return propellor_design(V, df['Pe'], df['L'], df['B'], df['T'], df['Cb'], df['Cp'], df['LCB'])
    
    hullsdf['Pd'] = hullsdf.apply(Pd, axis=1)
    
    eta_T = 0.98
    x = 0
    
    def Pi(df):
        try:
            return (df['Pd'] * (1 + x))/(eta_T)
        except:
            return 'NaN'
        
    hullsdf['Pi'] = hullsdf.apply(Pi, axis=1)
    
    def eedi(df):
        return EEDI(df['Pi'], DWT, V, fuel)
    
    
    hullsdf['EEDI'] = hullsdf.apply(eedi, axis=1)
    hullsdf = hullsdf.drop(hullsdf[hullsdf.Pd == 10**10].index)
    
    
    """-------------------------- Plots ------------------------------------"""
    """Ref EEDI"""
    
    
    min_idx = hullsdf['EEDI'].idxmin()
    max_idx = hullsdf['EEDI'].idxmax()
     
    L_min, B_min, T_min, Cb_min = hullsdf.at[min_idx, 'L'], hullsdf.at[min_idx, 'B'], hullsdf.at[min_idx, 'T'], hullsdf.at[min_idx, 'Cb']
    L_max, B_max, T_max, Cb_max = hullsdf.at[max_idx, 'L'], hullsdf.at[max_idx, 'B'], hullsdf.at[max_idx, 'T'], hullsdf.at[max_idx, 'Cb']
    
    EEDI_plot(hullsdf['EEDI'].min(), hullsdf['EEDI'].max(), DWT, vessel_type, L_min, L_max, B_min, B_max, T_min, T_max, Cb_min, Cb_max)
    

    hullsdf['L/D^(1/3)'] = hullsdf['L']/hullsdf['D^(1/3)']
    hullsdf['L/B'] = hullsdf['L']/hullsdf['B']
    hullsdf['B/T'] = hullsdf['B']/hullsdf['T']
    
    colors = ['lightgreen', 'gold', 'darkorange', 'red']
    labels = ["q1", "q2", "q3", "q4"]

    if dimensions != []:

        q1 = float(hullsdf['EEDI'].quantile([.25]))
        m = float(hullsdf['EEDI'].quantile([.5]))
        q3 = float(hullsdf['EEDI'].quantile([.75]))

        print(q1)

        q1_df = hullsdf[hullsdf['EEDI']<q1]
        q2_df = hullsdf.loc[(hullsdf['EEDI'] >= q1) & (hullsdf['EEDI'] <= m)]
        q3_df = hullsdf.loc[(hullsdf['EEDI'] >= m) & (hullsdf['EEDI'] <= q3)]
        q4_df = hullsdf[hullsdf['EEDI']>q3]

        m = 2*n 

        for i in dimensions:
            if i == "length":

                L1 = q1_df["L"].to_numpy()
                L2 = q2_df["L"].to_numpy()
                L3 = q3_df["L"].to_numpy()
                L4 = q4_df["L"].to_numpy()

                fig, ax = plt.subplots()

                ax.hist([L1, L2, L3, L4], m, histtype='bar', color=colors, stacked=True, label=labels)
                ax.legend()
                ax.set_xlabel('Length(m)')
                ax.set_ylabel('N Hulls')
                ax.set_title('Distribution of Length based on EEDI')
                
                fig.savefig(f'static/images/Stat_L {DWT} {V} {n}.png')

            elif i == "breadth":

                B1 = q1_df["B"].to_numpy()
                B2 = q2_df["B"].to_numpy()
                B3 = q3_df["B"].to_numpy()
                B4 = q4_df["B"].to_numpy()

                fig, ax = plt.subplots()

                ax.hist([B1, B2, B3, B4], m, histtype='bar', color=colors, stacked=True, label=labels)
                ax.legend()
                ax.set_xlabel('Breadth(m)')
                ax.set_ylabel('N Hulls')
                ax.set_title('Distribution of Breadth based on EEDI')
                
                fig.savefig(f'static/images/Stat_B {DWT} {V} {n}.png')

            elif i == "draught":

                T1 = q1_df["T"].to_numpy()
                T2 = q2_df["T"].to_numpy()
                T3 = q3_df["T"].to_numpy()
                T4 = q4_df["T"].to_numpy()

                fig, ax = plt.subplots()

                ax.hist([T1, T2, T3, T4], m, histtype='bar', color=colors, stacked=True, label=labels)
                ax.legend()
                ax.set_xlabel('Draught(m)')
                ax.set_ylabel('N Hulls')
                ax.set_title('Distribution of Draught based on EEDI')
                
                fig.savefig(f'static/images/Stat_T {DWT} {V} {n}.png')

            elif i == "cb":

                Cb1 = q1_df["Cb"].to_numpy()
                Cb2 = q2_df["Cb"].to_numpy()
                Cb3 = q3_df["Cb"].to_numpy()
                Cb4 = q4_df["Cb"].to_numpy()

                fig, ax = plt.subplots()

                ax.hist([Cb1, Cb2, Cb3, Cb4], n, histtype='bar', color=colors, stacked=True, label=labels)
                ax.legend()
                ax.set_xlabel('Block Coefficient')
                ax.set_ylabel('N Hulls')
                ax.set_title('Distribution of Block Coefficient based on EEDI')
                
                fig.savefig(f'static/images/Stat_Cb {DWT} {V} {n}.png')

            elif i == "ld":

                ld1 = q1_df["L/D^(1/3)"].to_numpy()
                ld2 = q2_df["L/D^(1/3)"].to_numpy()
                ld3 = q3_df["L/D^(1/3)"].to_numpy()
                ld4 = q4_df["L/D^(1/3)"].to_numpy()

                fig, ax = plt.subplots()

                ax.hist([ld1, ld2, ld3, ld4], n, histtype='bar', color=colors, stacked=True, label=labels)
                ax.legend()
                ax.set_xlabel('L/D^(1/3)')
                ax.set_ylabel('N Hulls')
                ax.set_title('Distribution of Slenderness based on EEDI')
                
                fig.savefig(f'static/images/Stat_LD {DWT} {V} {n}.png')

            elif i == "lb":

                lb1 = q1_df["L/B"].to_numpy()
                lb2 = q2_df["L/B"].to_numpy()
                lb3 = q3_df["L/B"].to_numpy()
                lb4 = q4_df["L/B"].to_numpy()

                fig, ax = plt.subplots()

                ax.hist([lb1, lb2, lb3, lb4], n, histtype='bar', color=colors, stacked=True, label=labels)
                ax.legend()
                ax.set_xlabel('Length/Breadth')
                ax.set_ylabel('N Hulls')
                ax.set_title('Repartition of L/B based on EEDI')
                
                fig.savefig(f'static/images/Stat_LB {DWT} {V} {n}.png')

            elif i == "bt":

                bt1 = q1_df["B/T"].to_numpy()
                bt2 = q2_df["B/T"].to_numpy()
                bt3 = q3_df["B/T"].to_numpy()
                bt4 = q4_df["B/T"].to_numpy()

                fig, ax = plt.subplots()

                ax.hist([bt1, bt2, bt3, bt4], n, histtype='bar', color=colors, stacked=True, label=labels)
                ax.legend()
                ax.set_xlabel('B/T')
                ax.set_ylabel('N Hulls')
                ax.set_title('Repartition of B/T based on EEDI')
                
                fig.savefig(f'static/images/Stat_BT {DWT} {V} {n}.png')



    
    return hullsdf[['L', 'B','T','Cb', 'EEDI', 'Pi']]

def eedi_study(V, DWT, vessel_type, fuel, years):
    
    hulls = shipdesign(V, DWT, vessel_type, fuel)
    eedi = time_eedi(DWT, years)
    
    percentage_eedi(eedi, years, hulls, DWT)
    
    hulls1 = []
    
    L0, B0, T0, Cb0 = opt_hull_eedi(eedi[0], hulls)
    below_eedi_hulls = hulls[hulls['EEDI']<eedi[0]]
    size = below_eedi_hulls.size
    
    for i in range(len(eedi)):
        try:
            L, B, T, Cb = opt_hull_eedi(eedi[i], hulls)
            hulls1.append([years[i], eedi[i], round(L, 1), round(B, 2), round(T, 2), round(Cb, 2) ])
            Lbelow = below_eedi_hulls[below_eedi_hulls['L']<L].size 
            perc_L = Lbelow/size
            print(perc_L)
        except:
            pass
            #hulls1.append([years[i], eedi[i], "NaN", "NaN", "NaN", "NaN", "NaN"])
        
    hulls2 = pd.DataFrame(hulls1, columns=['Year' , 'Actual EEDI','L', 'B',  'T', 'Cb' ])

    return hulls2

def air_lube_est(V, DWT, vessel_type, fuel, cf, area):
    
    data = []
    
    air_lube_cf = [0, cf]
    
    for i in air_lube_cf:
        hulls = shipdesign(V, DWT, vessel_type, fuel, {"Cf_lube":i, "airlubearea": area})
        
        Pmin = hulls.Pi.min()
        Pmax = hulls.Pi.max()
        Lmin = hulls.L.min()
        Lmax = hulls.L.max()
        
        
        q1 = float(hulls['Pi'].quantile([.25]))
        m = float(hulls['Pi'].quantile([.5]))
        q3 = float(hulls['Pi'].quantile([.75]))
        
        data.append([i, Pmin, q1, m, q3, Pmax, float(hulls.shape[0])])
        
        
        q1_df = hulls[hulls['Pi']<q1]
        q2_df = hulls.loc[(hulls['Pi'] >= q1) & (hulls['Pi'] <= m)]
        q3_df = hulls.loc[(hulls['Pi'] >= m) & (hulls['Pi'] <= q3)]
        q4_df = hulls[hulls['Pi']>q3]
        
        L1 = q1_df["L"].to_numpy()
        L2 = q2_df["L"].to_numpy()
        L3 = q3_df["L"].to_numpy()
        L4 = q4_df["L"].to_numpy()
        
        B1 = q1_df["B"].to_numpy()
        B2 = q2_df["B"].to_numpy()
        B3 = q3_df["B"].to_numpy()
        B4 = q4_df["B"].to_numpy()
        
        T1 = q1_df["T"].to_numpy()
        T2 = q2_df["T"].to_numpy()
        T3 = q3_df["T"].to_numpy()
        T4 = q4_df["T"].to_numpy()
        
        Cb1 = q1_df["Cb"].to_numpy()
        Cb2 = q2_df["Cb"].to_numpy()
        Cb3 = q3_df["Cb"].to_numpy()
        Cb4 = q4_df["Cb"].to_numpy()
        
        L_data = [L1, L2, L3, L4]
        B_data = [B1, B2, B3, B4]
        T_data = [T1, T2, T3, T4]
        Cb_data = [Cb1, Cb2, Cb3, Cb4]
        
        
        colors = ['lightgreen', 'gold', 'darkorange', 'red']
        labels = ["q1", "q2", "q3", "q4"]
        
        fig, ax = plt.subplots()
        
        ax.hist(L_data, 20, histtype='bar', color=colors,stacked=True, label=labels)
        ax.legend()
        ax.set_xlabel('L(m)')
        ax.set_ylabel('N Hulls')
        
        fig.savefig(f"static/images/ALS/{V} {DWT} Cf={i} A={area} L.png")
        
        fig, ax = plt.subplots()
        
        ax.hist(B_data, 20, histtype='bar', color=colors,stacked=True, label=labels)
        ax.legend()
        ax.set_xlabel('B(m)')
        ax.set_ylabel('N Hulls')
        
        fig.savefig(f"static/images/ALS/{V} {DWT} Cf={i} A={area} B.png")
        
        fig, ax = plt.subplots()
        
        ax.hist(T_data, 20, histtype='bar', color=colors,stacked=True, label=labels)
        ax.legend()
        ax.set_xlabel('T(m)')
        ax.set_ylabel('N Hulls')
        
        fig.savefig(f"static/images/ALS/{V} {DWT} Cf={i} A={area} T.png")
        
        fig, ax = plt.subplots()
        
        ax.hist(Cb_data, 20, histtype='bar', color=colors,stacked=True, label=labels)
        ax.legend()
        ax.set_xlabel('Cb')
        ax.set_ylabel('N Hulls')
        
        fig.savefig(f"static/images/ALS/{V} {DWT} Cf={i} A={area} Cb.png")
        
    data1 = pd.DataFrame(data, columns = ['Cf reduction', 'Minimum Pi', '1st Quartile', 'Median', '2nd Quartile', 'Maximum Pi', 'N'])
   
    return data1


