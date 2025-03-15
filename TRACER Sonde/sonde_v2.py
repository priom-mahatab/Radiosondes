#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 13:14:41 2025

@author: zmahatab
"""
import os
from metpy.units import units
import matplotlib.pyplot as plt
import metpy.plots as plots
import numpy as np
import metpy.calc as mpcalc
import pandas as pd
import matplotlib.gridspec as gridspec
from metpy.plots import Hodograph

import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('png', 'jpeg')

folder_path = './ozonesonde data files/'
output_folder = './Plots/'

txt_files = [f for f in os.listdir(folder_path) if f.endswith(".txt")]

for file in txt_files:
    print("Processing file", file, "\n")
    file_path = folder_path + file
    
    launch_date = ""
    launch_time_utc = ""
    launch_time_short = ""
    
    with open(file_path, 'r') as file:
        
        first_line = file.readline().strip()  
        if first_line[:2].isdigit():  
            n_header_lines = int(first_line[:2])
            
        for line in file:
            line = line.strip()
            
            parts = line.split(",")
            
            if parts[0] == "2022":
                year = parts[0].strip()
                month = parts[1].strip()
                day = parts[2].strip()
                launch_date = year+month+day
                launch_date_for_title = f"{month}-{day}-{year}"
                continue
            
            if line.startswith("DATA_INFO") and "=" in line:
                launch_time_utc = line.split("=")[-1].strip()  # Get the time after '='
                launch_time_short = launch_time_utc[:2] + launch_time_utc[3:5]
                break
                
            
                        
    
    df = pd.read_csv(file_path, delimiter=",", skiprows=n_header_lines-1)
    
    # remove leading and trailing whitespace from column names
    df.columns = df.columns.str.strip()
    
    print(df.head(), "\n")
    print(df.columns, "\n")
    
    # np.where(condition, value_if_true, value_if_false)
    df['WindSpeed_mps'] = np.where(df['WindSpeed_mps'] == -99999.0, np.nan, df['WindSpeed_mps'])
    df['WindDirection_deg'] = np.where(df['WindDirection_deg'] == -99999.0, np.nan, df['WindDirection_deg'])
    
    Pressure_hPa_list = []
    Temp_degC_list = []
    Frostpoint_degC_list = []
    WindSpeed_mps_list = []
    WindDirection_deg_list = []
    Altitude_km_list = []
    
    # Create lists that do not have any values that are flagged as -99999 (missing data)
    for i in range(len(df['Pressure_hPa'])):
        if df['Pressure_hPa'][i] < -9000:
            continue
        elif df['Temp_degC'][i] < -9000:
            continue
        # elif df['Frostpoint_degC'][i] < -9000:
        #     continue
        # elif df['WindSpeed_mps'][i] < -9000:
        #      df['WindSpeed_mps'][i] = np.nan
        #     continue
        # elif df['WindDirection_deg'][i] < -9000:
        #      df['WindDirection_deg'][i] = np.nan
        #     continue
        elif df['Altitude_km'][i] < -9000:
            continue
        
        Pressure_hPa_list.append(df['Pressure_hPa'][i])
        Temp_degC_list.append(df['Temp_degC'][i])
        Frostpoint_degC_list.append(df['Frostpoint_degC'][i])
        WindSpeed_mps_list.append(df['WindSpeed_mps'][i])
        WindDirection_deg_list.append(df['WindDirection_deg'][i])
        Altitude_km_list.append(df['Altitude_km'][i])
        
    Pressure_hPa_array = np.asarray(Pressure_hPa_list)
    Temp_degC_array = np.asarray(Temp_degC_list)
    Frostpoint_degC_array = np.asarray(Frostpoint_degC_list)
    WindSpeed_mps_array = np.asarray(WindSpeed_mps_list)
    WindDirection_deg_array = np.asarray(WindDirection_deg_list)
    Altitude_km_array = np.asarray(Altitude_km_list)
    
    # Determine when the launch occured and when goes past 100 hPa
    for i in range(len(Pressure_hPa_array)):
        pressure = Pressure_hPa_array[i]
        if i == 0:
            pressure_old = pressure
        if abs(pressure - pressure_old) < 0.05 and i < 100:
            launch_cell = i
            if launch_cell == 100:
                print("Check launch cell in file", file)
    
        if pressure < 90:
            end_cell = i
            break
    print("The launch cell is", launch_cell, " and the end cell is", end_cell, "\n")
    
    Pressure_hPa_array = Pressure_hPa_array[launch_cell:end_cell]
    Temp_degC_array = Temp_degC_array[launch_cell:end_cell]
    Frostpoint_degC_array = Frostpoint_degC_array[launch_cell:end_cell]
    WindSpeed_mps_array = WindSpeed_mps_array[launch_cell:end_cell]
    WindDirection_deg_array = WindDirection_deg_array[launch_cell:end_cell]
    Altitude_km_array = Altitude_km_array[launch_cell:end_cell]
    
    p = Pressure_hPa_array * units.hPa
    T = Temp_degC_array * units.degC
    Td = Frostpoint_degC_array * units.degC
    WindSpeed= WindSpeed_mps_array * units.mps
    WindDir = WindDirection_deg_array * units.deg
    heights = Altitude_km_array * units.km
    
    # Determine the u and v components of the wind
    u, v = mpcalc.wind_components(WindSpeed, WindDir)
    
    try:
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
        el_pressure, el_temperature = mpcalc.el(p, T, Td)
        lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td)
        surface_cape, surface_cin = mpcalc.surface_based_cape_cin(p, T, Td)
        mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(p, T, Td)
        parcel_path = mpcalc.parcel_profile(p, T[0], Td[0])
    except:
        print("Skipping file:", file, "due to error")
        continue
    
    print("LCL pressure", lcl_pressure)
    print("EL pressure", el_pressure)
    print("LFC pressure", lfc_pressure)
    print("Surface CAPE/CIN", surface_cape, surface_cin)
    print("Most unstable CAPE/CIN", mu_cape, mu_cin)
    
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[5,2.5])
    gs_bottom = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=gs[1], width_ratios=[3,2])
    
    
    # SkewT
    if launch_time_utc[:2] == "10" or launch_time_utc[:2] == "11":
        
        skew = plots.SkewT(fig, subplot=gs[0])
        skew.ax.set_xlim(-90,30)
        skew.ax.set_ylim(1050,100)
        skew.plot(p, T, 'red', label="Temperature")
        skew.plot(p, Td, 'green', label="Dew Point")
        
        # hodograph
        ax_hod = fig.add_subplot(gs_bottom[0])
        hod = Hodograph(ax_hod, component_range=40)
        hod.add_grid(increment=10)
        hod.add_grid(increment=20, color='tab:orange', linestyle='-')
        hod.plot_colormapped(u, v, WindSpeed, linewidth=2)
        
        interval = np.logspace(2, 3) * units.hPa
        idx = mpcalc.resample_nn_1d(p, interval)
        skew.plot_barbs(p[idx][::2], u[idx][::2], v[idx][::2])
        
        skew.plot_dry_adiabats()
        skew.plot_moist_adiabats()
        skew.plot_mixing_lines()
        
        skew.plot(p, parcel_path, color='k')
        skew.shade_cape(p, T, parcel_path)
        skew.shade_cin(p, T, parcel_path)
        
        if lcl_pressure:
            skew.ax.axhline(lcl_pressure, color ='orange', label='LCL Pressure')
        if lfc_pressure:
            skew.ax.axhline(lfc_pressure, color='purple', label='LFC Pressure')
        if el_pressure:
            skew.ax.axhline(el_pressure, color='blue', label='Equilibrium level')
        
    else:
        skew = plots.SkewT(fig, subplot=gs[0])
        skew.ax.set_xlim(-90,30)
        skew.ax.set_ylim(1050,100)
        skew.plot(p, T, 'red', label="Temperature", linestyle="dashed", alpha=0.8)
        skew.plot(p, Td, 'green', label="Dew Point", linestyle="dashed", alpha=0.8)
        
        # hodograph
        ax_hod = fig.add_subplot(gs_bottom[0])
        hod = Hodograph(ax_hod, component_range=40)
        hod.add_grid(increment=10)
        hod.add_grid(increment=20, color='tab:orange', linestyle='-')
        hod.plot_colormapped(u, v, WindSpeed, linewidth=2)
        
        interval = np.logspace(2, 3) * units.hPa
        idx = mpcalc.resample_nn_1d(p, interval)
        skew.plot_barbs(p[idx][::2], u[idx][::2], v[idx][::2])
        
        skew.plot_dry_adiabats(alpha=0.5)
        skew.plot_moist_adiabats(alpha=0.5)
        skew.plot_mixing_lines(alpha=0.5)
        
        skew.plot(p, parcel_path, color='k', linestyle="dashed", alpha=0.8)
        skew.shade_cape(p, T, parcel_path)
        skew.shade_cin(p, T, parcel_path)
        
        if lcl_pressure:
            skew.ax.axhline(lcl_pressure, color ='orange', label='LCL Pressure', linestyle="dashed", alpha=0.8)
        if lfc_pressure:
            skew.ax.axhline(lfc_pressure, color='purple', label='LFC Pressure', linestyle="dashed", alpha=0.8)
        if el_pressure:
            skew.ax.axhline(el_pressure, color='blue', label='Equilibrium level', linestyle="dashed", alpha=0.8)
    
    
    
    
    # plots.add_timestamp(skew.ax, y=1.05, fontsize=18)
    skew.ax.legend(loc='upper right', bbox_to_anchor=(0.95, 1), fontsize=10, frameon=True)
    
    # stats axis
    ax_stats = fig.add_subplot(gs_bottom[1])
    ax_stats.set_xticks([])
    ax_stats.set_yticks([])
    ax_stats.spines['top'].set_visible(False)
    ax_stats.spines['bottom'].set_visible(False)
    ax_stats.spines['left'].set_visible(False)
    ax_stats.spines['right'].set_visible(False)
    
    # stats dictionary
    stats_values = {
        r"MUCAPE (J$\,$kg$^{-1}$)": mu_cape.magnitude,
        r"MUCIN (J$\,$kg$^{-1}$)": mu_cin.magnitude,
        r"Surface CAPE (J$\,$kg$^{-1}$)": surface_cape.magnitude,
        r"Surface CIN (J$\,$kg$^{-1}$)": surface_cin.magnitude,
        "LCL Pressure (hPa)": lcl_pressure.magnitude,
        r"LCL Temperature ($^\circ$C)": lcl_temperature.magnitude,
        "LFC Pressure (hPa)": lfc_pressure.magnitude,
        r"LFC Temperature ($^\circ$C)": lfc_temperature.magnitude,
        "Equilibrium Level (hPa)": el_temperature.magnitude,
        r"Equilibrium Level ($^\circ$C)": el_pressure.magnitude
    }
    
    # value location
    x_left = -0.3
    x_right = 0.4
    
    y_start = 0.95
    y_spacing = 0.1
    
    for i, (label, value) in enumerate(stats_values.items()):
        ax_stats.text(x_left, y_start - i * y_spacing, f"{label}:", fontsize=10, ha='left', va='top', family='monospace') # labels
        ax_stats.text(x_right, y_start - i * y_spacing, f"{value:.2f}", fontsize=10, ha='right', va='top',family='monospace') # values
    
    plt.tight_layout()
    
    ax_hod.set_position([0.19, 0.05, 0.25, 0.25])

    skew.ax.set_title('TRACER Sonde ' + launch_date_for_title + " " + launch_time_short[0:2] + ":" + launch_time_short[2:], fontsize=10)
    fig_filename = output_folder + "TRACERSonde." + launch_date + "." + launch_time_short + '.stats.png'
          
    plt.show()
    
    fig.savefig(fig_filename, bbox_inches='tight', dpi=300)
