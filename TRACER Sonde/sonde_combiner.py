#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 2:15:41 2025

@author: zmahatab
"""

import os
from metpy.units import units
import matplotlib.pyplot as plt
import metpy.plots as plots
import numpy as np
import metpy.calc as mpcalc
import pandas as pd

folder_path = '/Users/zmahatab/Desktop/MetPy/TRACER Sonde/ozonesonde data files/'
output_folder = '/Users/zmahatab/Desktop/MetPy/TRACER Sonde/Plots/Combined'

def get_launch_info(file_path):
    launch_date = ""
    launch_time_utc = ""
    launch_time_short = ""
    
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        if first_line[:2].isdigit():
            n_header_lines = int(first_line[:2])
        else:
            return None, None, None, None
        
        for line in file:
            line = line.strip()
            parts = line.split(",")
            if parts[0] == "2022":
                year, month, day = parts[:3]
                launch_date = year.strip() + month.strip() + day.strip()
                launch_date_for_title = f"{month}-{day}-{year}"
                
            if line.startswith("DATA_INFO") and "=" in line:
                launch_time_utc = line.split("=")[-1].strip()
                launch_time_short = launch_time_utc[:2] + launch_time_utc[3:5]
                break
    
    return launch_date, launch_time_utc, launch_time_short, n_header_lines

def check_pressure(p):
    pressure_increases = False
    for i in range(len(p)-1):
        if p[i+1] >= p[i]:
            print("Pressure increases", i, p[i], i+1, p[i+1])
            pressure_increases = True
    return pressure_increases


def process_file(file_path, n_header_lines):
    df = pd.read_csv(file_path, delimiter=",", skiprows=n_header_lines - 1)
    df.columns = df.columns.str.strip()
    
    for col in ['WindSpeed_mps', 'WindDirection_deg']:
        if col in df.columns:
            df[col] = df[col].replace(-99999.0, np.nan)
    
    valid_rows = (df['Pressure_hPa'] > -9000) & (df['Temp_degC'] > -9000) & (df['Altitude_km'] > -9000)
    df = df[valid_rows]
    
    return df

def create_skewT(launch_date_for_title, df1, df2, output_path):
    fig = plt.figure(figsize=(12, 8))
    skew = plots.SkewT(fig)
    skew.ax.set_xlim(-90, 30)
    skew.ax.set_ylim(1050, 100)
    
    # First file data
    p = df1['Pressure_hPa'].values * units.hPa
    T = df1['Temp_degC'].values * units.degC
    Td = df1['Frostpoint_degC'].values * units.degC
    u, v = mpcalc.wind_components(df1['WindSpeed_mps'].values * units.mps, 
                                  df1['WindDirection_deg'].values * units.deg)
    
    coarsening = 0
    pressure_increases_check = check_pressure(p)
    while pressure_increases_check:
        p = p[::2]
        T = T[::2]
        Td = Td[::2]
        u = u[::2]
        v = v[::2]
        coarsening += 1
        pressure_increases_check = check_pressure(p)
    print("Exiting pressure check")
        
    try:
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
        el_pressure, el_temperature = mpcalc.el(p, T, Td)
        lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td)
        surface_cape, surface_cin = mpcalc.surface_based_cape_cin(p, T, Td)
        mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(p, T, Td)
        parcel_path = mpcalc.parcel_profile(p, T[0], Td[0])
    except:
        print("Skipping file: due to error")
        
    if lcl_pressure:
        skew.ax.axhline(lcl_pressure, color ='orange', label='LCL Pressure')
    if lfc_pressure:
        skew.ax.axhline(lfc_pressure, color='purple', label='LFC Pressure')
    if el_pressure:
        skew.ax.axhline(el_pressure, color='blue', label='Equilibrium level')
    
    # Second file data
    p2 = df2['Pressure_hPa'].values * units.hPa
    T2 = df2['Temp_degC'].values * units.degC
    Td2 = df2['Frostpoint_degC'].values * units.degC
    u2, v2 = mpcalc.wind_components(df2['WindSpeed_mps'].values * units.mps, 
                                    df2['WindDirection_deg'].values * units.deg)
    
    coarsening2 = 0
    pressure2_increases_check = check_pressure(p2)
    
    
    while pressure2_increases_check:
        p2 = p2[::2]
        T2 = T2[::2]
        Td2 = Td2[::2]
        u2 = u2[::2]
        v2 = v2[::2]
        coarsening2 += 1
        pressure2_increases_check = check_pressure(p2)
    
    print("Exiting pressure2 check")
        
    try:
        lcl_pressure2, lcl_temperature2 = mpcalc.lcl(p2[0], T2[0], Td2[0])
        el_pressure2, el_temperature2 = mpcalc.el(p2, T2, Td2)
        lfc_pressure2, lfc_temperature2 = mpcalc.lfc(p2, T2, Td2)
        surface_cape2, surface_cin2 = mpcalc.surface_based_cape_cin(p2, T2, Td2)
        mu_cape2, mu_cin2 = mpcalc.most_unstable_cape_cin(p2, T2, Td2)
        parcel_path2 = mpcalc.parcel_profile(p2, T2[0], Td2[0])
    except:
        print("Skipping file: due to error")
    
    if lcl_pressure2:
        skew.ax.axhline(lcl_pressure2, color ='orange', label='LCL Pressure', linestyle="dashed", alpha=0.6)
    if lfc_pressure2:
        skew.ax.axhline(lfc_pressure2, color='purple', label='LFC Pressure', linestyle="dashed", alpha=0.6)
    if el_pressure2:
        skew.ax.axhline(el_pressure2, color='blue', label='Equilibrium level', linestyle="dashed", alpha=0.6)
    
    # Plot first profile
    skew.plot(p, T, 'red', label="Temperature 1")
    skew.plot(p, Td, 'green', label="Dew Point 1")
    
    skew.plot(p, parcel_path, color='k')
    skew.plot(p2, parcel_path2, color='k', linestyle="--", alpha=0.6)
    skew.shade_cape(p, T, parcel_path)
    skew.shade_cin(p, T, parcel_path)
    
    interval = np.logspace(2, 3) * units.hPa
    idx = mpcalc.resample_nn_1d(p, interval)
    idx2 = mpcalc.resample_nn_1d(p2, interval)
    
    skew.plot_barbs(p[idx][::2], u[idx][::2], v[idx][::2])
    
    # Plot second profile with different opacity
    skew.plot(p2, T2, 'red', linestyle="--", alpha=0.6, label="Temperature 2")
    skew.plot(p2, Td2, 'green', linestyle="--", alpha=0.6, label="Dew Point 2")
    skew.plot_barbs(p2[idx2][::2], u2[idx2][::2], v2[idx2][::2], color="red")
    
    # Add dry adiabats, moist adiabats, and mixing lines
    # skew.plot_dry_adiabats()
    # skew.plot_moist_adiabats()
    # skew.plot_mixing_lines()

    # Final formatting
    skew.ax.set_title(f'TRACER Sonde {launch_date_for_title}', fontsize=12)
    #skew.ax.legend(loc='upper right',bbox_to_anchor=(0.95, 1), fontsize=10, frameon=True)
    
    # Save figure
    fig.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)

def combined_skewT(folder_path, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    txt_files = [f for f in os.listdir(folder_path) if f.endswith(".txt")]
    
    for file1 in txt_files:
        file_path1 = os.path.join(folder_path, file1)
        launch_date1, launch_time1, time_short1, n_header_lines1 = get_launch_info(file_path1)
        
        if not launch_date1:
            print("Launch date does not exist for {file1}")
            continue
        
        if launch_time1[:2] in ["10", "11"]:
            
            print("Time found!")
            df1 = process_file(file_path1, n_header_lines1)
            
            for file2 in txt_files:
                if file1 == file2:
                    continue
                
                file_path2 = os.path.join(folder_path, file2)
                launch_date2, launch_time2, time_short2, n_header_lines2 = get_launch_info(file_path2)
                if not launch_date2:
                    print("Launch date does not exist for {file2}")
                    continue
                
                if launch_date1 == launch_date2 and launch_time1 != launch_time2:
                    df2 = process_file(file_path2, n_header_lines2)
                    output_path = os.path.join(output_folder, f"TRACERSonde.{launch_date1}.png")
                    
                    print(f"Creating plot for {file1} and {file2}")
                    
                    try:
                        create_skewT(f"{launch_date1[:4]}-{launch_date1[4:6]}-{launch_date1[6:]}", df1, df2, output_path)
                    except Exception as e:
                        print(f"Error creating plot for {file1} and {file2}: {e}")

# Run function
combined_skewT(folder_path, output_folder)
