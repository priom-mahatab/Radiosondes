import os
from datetime import datetime
from metpy.units import units
import matplotlib.gridspec as gridspec
import pandas as pd
import matplotlib.pyplot as plt
import metpy.plots as plots
from metpy.plots import SkewT, Hodograph
import metpy.calc as mpcalc
from metpy.calc import lcl, parcel_profile, cape_cin, lifted_index
import numpy as np
%matplotlib inline

input_folder = "/Users/zmahatab/Desktop/MetPy/TRACER_M1/Data"
output_folder = "/Users/zmahatab/Desktop/MetPy/TRACER_M1/Plots"

os.makedirs(output_folder, exist_ok=True)

def skewT_generator_11(input_folder, output_folder):
    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            parts = filename.split(".")
            if len(parts) >= 5: # ensuring the file has correct naming format
                date_part = parts[2] #YYYYMMDD
                time_part = parts[3] # HHMMSS
                hour_part = parts[3][:2] # HH
    
                if hour_part == "11":
                    file_path = os.path.join(input_folder, filename)
    
                    # Format date: YYYYMMDD -> MM-DD-YYYY
                    formatted_date = f"{date_part[4:6]}-{date_part[6:8]}-{date_part[:4]}"
    
                    # Format time: HHMMSS -> HH:MM
                    formatted_time = f"{time_part[:2]}:{time_part[2:4]}"
                    print(f"Processing {filename}")
    
                    data = pd.read_csv(file_path)
    
                    # extracting data
                    pressure = data['pres'].values * units.hPa
                    dry_temp = data['tdry'].values * units.degC
                    dew_point = data['dp'].values * units.degC
                    wind_speed = data['wspd'].values * units.knot
                    u_wind = data['u_wind'].values * units.knot
                    v_wind = data['v_wind'].values * units.knot
                    height = data['alt'].values * units.metre
    
                    if pressure.size == 0 or dry_temp.size == 0 or dew_point.size == 0:
                        print(f"Skipping {filename} due to empty data.")
                        continue
    
                    # initializing figure and slots
                    fig = plt.figure(figsize=(12,8))
                    gs = gridspec.GridSpec(2, 1, height_ratios=[5,2.5])
                    gs_bottom = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=gs[1], width_ratios=[3,2])
    
                    # SkewT
                    skew = plots.SkewT(fig, subplot=gs[0])
                    skew.ax.set_xlim(-90,30)
                    skew.ax.set_ylim(1050,100)
    
                    # hodograph
                    ax_hod = fig.add_subplot(gs_bottom[0])
                    hod = Hodograph(ax_hod, component_range=40)
                    hod.add_grid(increment=10)
                    hod.add_grid(increment=20, color='tab:orange', linestyle='-')
    
                    # plotting
                    skew.plot(pressure, dry_temp, 'red', label='Temperature')
                    skew.plot(pressure, dew_point, 'green', label='Dew Point')
    
                    hod.plot_colormapped(u_wind, v_wind, wind_speed, linewidth=2)
    
                    # fixing wind barbs
                    interval = np.logspace(2, 3) * units.hPa
                    idx = mpcalc.resample_nn_1d(pressure, interval)
                    skew.plot_barbs(pressure[idx][::2], u_wind[idx][::2], v_wind[idx][::2])
                    
                    # parcel path
                    parcel_path = mpcalc.parcel_profile(pressure, dry_temp[0], dew_point[0])
           
                    # stats
                    lfc_pressure, lfc_temp = mpcalc.lfc(pressure, dry_temp, dew_point)
                    lcl_pressure, lcl_temp = mpcalc.lcl(pressure[0], dry_temp[0], dew_point[0])
                    surface_cape, surface_cin = mpcalc.surface_based_cape_cin(pressure, dry_temp, dew_point)
                    li = lifted_index(pressure, dry_temp, parcel_path)
                    el = mpcalc.el(pressure, dry_temp, dew_point)
                    try:
                        mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(pressure, dry_temp, dew_point)
                    except ValueError as e:
                        print(f"Skipping file due to error: {e}")
                        continue
    
                    # fiducial lines
                    skew.plot_dry_adiabats()
                    skew.plot_moist_adiabats()
                    skew.plot_mixing_lines()
    
                    skew.plot(pressure, parcel_path, color='black')
                    skew.shade_cape(pressure, dry_temp, parcel_path)
                    skew.shade_cin(pressure, dry_temp, parcel_path)
                    
                    # horizontal lines
                    skew.ax.axhline(el[0], color="blue", label="Equilibrium level")
                    skew.ax.axhline(lfc_pressure, color="purple", label="LFC Pressure")
                    skew.ax.axhline(lcl_pressure, color="orange", label="LFC Pressure")
    
                    title = f"TRACER M1 {formatted_date} {formatted_time}"
                    skew.ax.set_title(title, fontsize=10)
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
                        r"LCL Temperature ($^\circ$C)": lcl_temp.magnitude,
                        "LFC Pressure (hPa)": lfc_pressure.magnitude,
                        r"LFC Temperature ($^\circ$C)": lfc_temp.magnitude,
                        "Equilibrium Level (hPa)": el[0].magnitude,
                        r"Equilibrium Level ($^\circ$C)": li.magnitude.item()
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
    
                    save_path = f"{output_folder}/TRACERM1.{date_part}.{time_part}.stats.png"
                    fig.savefig(save_path, dpi=300)


def skewT_generator_17(input_folder, output_folder):
    for filename in os.listdir(input_folder):
        if filename.endswith("csv"):
            parts = filename.split(".")
            if len(parts) >= 5: # ensuring the file has correct naming format
                date_part = parts[2] #YYYYMMDD
                time_part = parts[3] # HHMMSS
                hour_part = parts[3][:2] # HH
    
                if hour_part == "17":
                    file_path = os.path.join(input_folder, filename)

                    # Format date: YYYYMMDD -> MM-DD-YYYY
                    formatted_date = f"{date_part[4:6]}-{date_part[6:8]}-{date_part[:4]}"
    
                    # Format time: HHMMSS -> HH:MM
                    formatted_time = f"{time_part[:2]}:{time_part[2:4]}"
                    print(f"Processing {filename}")
    
                    data = pd.read_csv(file_path)
    
                    # extracting data
                    pressure = data['pres'].values * units.hPa
                    dry_temp = data['tdry'].values * units.degC
                    dew_point = data['dp'].values * units.degC
                    wind_speed = data['wspd'].values * units.knot
                    u_wind = data['u_wind'].values * units.knot
                    v_wind = data['v_wind'].values * units.knot
                    height = data['alt'].values * units.metre
    
                    if pressure.size == 0 or dry_temp.size == 0 or dew_point.size == 0:
                        print(f"Skipping {filename} due to empty data.")
                        continue

                    # initializing figure and slots
                    fig = plt.figure(figsize=(12,8))
                    gs = gridspec.GridSpec(2, 1, height_ratios=[5,2.5])
                    gs_bottom = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=gs[1], width_ratios=[3,2])
    
                    # SkewT
                    skew = plots.SkewT(fig, subplot=gs[0])
                    skew.ax.set_xlim(-90,30)
                    skew.ax.set_ylim(1050,100)
    
                    # hodograph
                    ax_hod = fig.add_subplot(gs_bottom[0])
                    hod = Hodograph(ax_hod, component_range=40)
                    hod.add_grid(increment=10)
                    hod.add_grid(increment=20, color='tab:orange', linestyle='-')
    
                    # plotting
                    skew.plot(pressure, dry_temp, 'red', label='Temperature', linestyle="dashed", alpha=0.6)
                    skew.plot(pressure, dew_point, 'green', label='Dew Point', linestyle="dashed", alpha=0.6)

                    hod.plot_colormapped(u_wind, v_wind, wind_speed, linewidth=2)
    
                    # fixing wind barbs
                    interval = np.logspace(2, 3) * units.hPa
                    idx = mpcalc.resample_nn_1d(pressure, interval)
                    skew.plot_barbs(pressure[idx][::2], u_wind[idx][::2], v_wind[idx][::2], color="red")
                    
                    # parcel path
                    parcel_path = mpcalc.parcel_profile(pressure, dry_temp[0], dew_point[0])
           
                    # stats
                    lfc_pressure, lfc_temp = mpcalc.lfc(pressure, dry_temp, dew_point)
                    lcl_pressure, lcl_temp = mpcalc.lcl(pressure[0], dry_temp[0], dew_point[0])
                    surface_cape, surface_cin = mpcalc.surface_based_cape_cin(pressure, dry_temp, dew_point)
                    li = lifted_index(pressure, dry_temp, parcel_path)
                    el = mpcalc.el(pressure, dry_temp, dew_point)
                    try:
                        mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(pressure, dry_temp, dew_point)
                    except ValueError as e:
                        print(f"Skipping file due to error: {e}")
                        continue

                    # fiducial lines
                    skew.plot_dry_adiabats(linestyle="dashed", alpha=0.6)
                    skew.plot_moist_adiabats(linestyle="dashed", alpha=0.6)
                    skew.plot_mixing_lines(linestyle="dashed", alpha=0.6)
    
                    skew.plot(pressure, parcel_path, color='black', linestyle="dashed", alpha=0.6)
                    skew.shade_cape(pressure, dry_temp, parcel_path)
                    skew.shade_cin(pressure, dry_temp, parcel_path)
                    
                    # horizontal lines
                    skew.ax.axhline(el[0], color="blue", label="Equilibrium level", linestyle="dashed", alpha=0.6)
                    skew.ax.axhline(lfc_pressure, color="purple", label="LFC Pressure", linestyle="dashed", alpha=0.6)
                    skew.ax.axhline(lcl_pressure, color="orange", label="LFC Pressure", linestyle="dashed", alpha=0.6)
    
                    title = f"TRACER M1 {formatted_date} {formatted_time}"
                    skew.ax.set_title(title, fontsize=10)
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
                        r"LCL Temperature ($^\circ$C)": lcl_temp.magnitude,
                        "LFC Pressure (hPa)": lfc_pressure.magnitude,
                        r"LFC Temperature ($^\circ$C)": lfc_temp.magnitude,
                        "Equilibrium Level (hPa)": el[0].magnitude,
                        r"Equilibrium Level ($^\circ$C)": li.magnitude.item()
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
    
                    save_path = f"{output_folder}/TRACERM1.{date_part}.{time_part}.stats.png"
                    fig.savefig(save_path, dpi=300)


def overlapping_skewTs(input_folder, output_folder):
    for file1 in os.listdir(input_folder):
        parts = filename.split(".")
        date_part = parts[2]
        time_part = parts[3]
        hour_part = parts[3][:2]

        if hour_part == "11":
            file_path1 = os.path.join(input_folder, file1)
            formatted_date = f"{date_part[4:6]}-{date_part[6:8]}-{date_part[:4]}"
            for file2 in os.listdir(input_folder):
                parts2 = file2.split(".")
                date_part2 = parts2[2]
                time_part2 = parts2[3][:2]
                
                if date_part == date_part2 and time_part2 == "17":
                    file_path2 = os.path.join(input_folder, file2)
                    print(f"Processing {file1} and {file2}")
                    data1 = pd.read_csv(file_path1)
                    data2 = pd.read_csv(file_path2)

                    # file1 data
                    pressure1 = data1['pres'].values * units.hPa
                    dry_temp1 = data1['tdry'].values * units.degC
                    dew_point1 = data1['dp'].values * units.degC
                    wind_speed1 = data1['wspd'].values * units.knot
                    u_wind1 = data1['u_wind'].values * units.knot
                    v_wind1 = data1['v_wind'].values * units.knot
                    height1 = data1['alt'].values * units.metre

                    # file2 data
                    pressure2 = data2['pres'].values * units.hPa
                    dry_temp2 = data2['tdry'].values * units.degC
                    dew_point2 = data2['dp'].values * units.degC
                    wind_speed2 = data2['wspd'].values * units.knot
                    u_wind2 = data2['u_wind'].values * units.knot
                    v_wind2 = data2['v_wind'].values * units.knot
                    height2 = data2['alt'].values * units.metre

                    # SkewT
                    fig = plt.figure(figsize=(12,8))
                    skew = plots.SkewT(fig)
                    skew.ax.set_xlim(-90,30)
                    skew.ax.set_ylim(1050,100)

                    skew.plot(pressure1, dry_temp1, 'red', label='Temperature1')
                    skew.plot(pressure1, dew_point1, 'green', label='Dew Point 1')
                    skew.plot(pressure2, dry_temp2, 'red', label='Temperature2', linestyle='dashed', alpha=0.6)
                    skew.plot(pressure2, dew_point2, 'green', label='Dew Point 2', linestyle='dashed', alpha=0.6)

                    #wind barbs
                    interval = np.logspace(2, 3) * units.hPa
                    idx1 = mpcalc.resample_nn_1d(pressure1, interval)
                    skew.plot_barbs(pressure1[idx1][::2], u_wind1[idx1][::2], v_wind1[idx1][::2])
                    idx2 = mpcalc.resample_nn_1d(pressure2, interval)
                    skew.plot_barbs(pressure2[idx2][::2], u_wind2[idx2][::2], v_wind2[idx2][::2], color="red")
                    
                    parcel_path1 = mpcalc.parcel_profile(pressure1, dry_temp1[0], dew_point1[0])
                    parcel_path2 = mpcalc.parcel_profile(pressure2, dry_temp2[0], dew_point2[0])

                    # stats
                    lfc_pressure1, lfc_temp1 = mpcalc.lfc(pressure1, dry_temp1, dew_point1)
                    lcl_pressure1, lcl_temp1 = mpcalc.lcl(pressure1[0], dry_temp1[0], dew_point1[0])
                    surface_cape1, surface_cin1 = mpcalc.surface_based_cape_cin(pressure1, dry_temp1, dew_point1)
                    li1 = lifted_index(pressure1, dry_temp1, parcel_path1)
                    el1 = mpcalc.el(pressure1, dry_temp1, dew_point1)
                    try:
                        mu_cape1, mu_cin1 = mpcalc.most_unstable_cape_cin(pressure1, dry_temp1, dew_point1)
                    except ValueError as e:
                        print(f"Skipping file due to error: {e}")
                        continue
                    
                    
                    lfc_pressure2, lfc_temp2 = mpcalc.lfc(pressure2, dry_temp2, dew_point2)
                    lcl_pressure2, lcl_temp2 = mpcalc.lcl(pressure2[0], dry_temp2[0], dew_point2[0])
                    surface_cape2, surface_cin2 = mpcalc.surface_based_cape_cin(pressure2, dry_temp2, dew_point2)
                    li2 = lifted_index(pressure2, dry_temp2, parcel_path2)
                    el2 = mpcalc.el(pressure2, dry_temp2, dew_point2)
                    try:
                        mu_cape2, mu_cin2 = mpcalc.most_unstable_cape_cin(pressure2, dry_temp2, dew_point2)
                    except ValueError as e:
                        print(f"Skipping file due to error: {e}")
                        continue

                    skew.plot(pressure1, parcel_path1, color='black')
                    skew.plot(pressure2, parcel_path2, color='black', alpha=0.6, linestyle='dashed')

                    # horizontal lines
                    skew.ax.axhline(el1[0], color='blue', label="Equilibrium level 1")
                    skew.ax.axhline(lfc_pressure1, color='purple', label="LFC Pressure 1")
                    skew.ax.axhline(lcl_pressure1, color='orange', label="LCL Pressure 1")
                    
                    skew.ax.axhline(el2[0], color='blue', label="Equilibrium level 2", linestyle='dashed', alpha=0.6)
                    skew.ax.axhline(lfc_pressure2, color='purple', label="LFC Pressure 2", linestyle='dashed', alpha=0.6)
                    skew.ax.axhline(lcl_pressure2, color='orange', label="LFC Pressure 2", linestyle='dashed', alpha=0.6)

                    skew.shade_cape(pressure1, dry_temp1, parcel_path1)
                    skew.shade_cin(pressure1, dry_temp1, parcel_path1)

                    title = f"TRACER M1 {formatted_date}"
                    skew.ax.set_title(title, fontsize=10)

                    save_path = f"{output_folder}/Combined/TRACERM1.{date_part}.png"
                    fig.savefig(save_path, dpi=300)