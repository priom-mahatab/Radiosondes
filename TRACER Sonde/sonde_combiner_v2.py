import os
from metpy.units import units
import matplotlib.pyplot as plt
import metpy.plots as plots
import numpy as np
import metpy.calc as mpcalc
from scipy.signal import medfilt
import pandas as pd

folder_path = '/Users/zmahatab/Desktop/MetPy/TRACER Sonde/ozonesonde data files/'
output_folder = '/Users/zmahatab/Desktop/MetPy/TRACER Sonde/Plots/Combined'

def combined_skewT(folder_path, output_folder):
    
    txt_files = [f for f in os.listdir(folder_path) if f.endswith(".txt")]
    
    for file in txt_files:
        file_path = folder_path + file
        
        launch_date = ""
        launch_time_utc = ""
        launch_time_short = ""
        
        with open(file_path, 'r') as file1:
            
            first_line = file1.readline().strip()  
            if first_line[:2].isdigit():  
                n_header_lines = int(first_line[:2])
                
            for line in file1:
                line = line.strip()
                
                parts = line.split(",")
                
                if parts[0] == "2022":
                    year = parts[0].strip()
                    month = parts[1].strip()
                    day = parts[2].strip()
                    launch_date = year+month+day
                    launch_date_for_title = f"{month}-{day}-{year}"
                    
                if line.startswith("DATA_INFO") and "=" in line:
                    launch_time_utc = line.split("=")[-1].strip()  # Get the time after '='
                    launch_time_short = launch_time_utc[:2] + launch_time_utc[3:5]
                    break

        for file2 in txt_files:
            file_path2 = folder_path + file2
            
            launch_date2 = ""
            launch_time_utc2 = ""
            launch_time_short2 = ""
            
            with open(file_path2, 'r') as file2:
                
                first_line2 = file2.readline().strip()  
                if first_line2[:2].isdigit():  
                    n_header_lines2 = int(first_line2[:2])
                    
                for line2 in file2:
                    line2 = line2.strip()
                    
                    parts2 = line2.split(",")
                    
                    if parts2[0] == "2022":
                        year2 = parts2[0].strip()
                        month2 = parts2[1].strip()
                        day2 = parts2[2].strip()
                        launch_date2 = year2+month2+day2
                        
                    if line.startswith("DATA_INFO") and "=" in line2:
                        launch_time_utc2 = line2.split("=")[-1].strip()  # Get the time after '='
                        launch_time_short2 = launch_time_utc2[:2] + launch_time_utc2[3:5] 
                        break
                    
            if launch_date == launch_date2: 
                print("Processing file", file1, " and", file2, "\n")
                df = pd.read_csv(file_path, delimiter=",", skiprows=n_header_lines-1)
                df2 = pd.read_csv(file_path2, delimiter=",", skiprows=n_header_lines2-1)
                
                # remove leading and trailing whitespace from column names
                df.columns = df.columns.str.strip()
                df2.columns = df2.columns.str.strip()
                
                df['WindSpeed_mps'] = np.where(df['WindSpeed_mps'] == -99999.0, np.nan, df['WindSpeed_mps'])
                df['WindDirection_deg'] = np.where(df['WindDirection_deg'] == -99999.0, np.nan, df['WindDirection_deg'])

                df2['WindSpeed_mps'] = np.where(df2['WindSpeed_mps'] == -99999.0, np.nan, df2['WindSpeed_mps'])
                df2['WindDirection_deg'] = np.where(df2['WindDirection_deg'] == -99999.0, np.nan, df2['WindDirection_deg'])
                
                # Apply median filter to smooth pressure data
                df['Pressure_hPa'] = medfilt(df['Pressure_hPa'], kernel_size=3)
                df2['Pressure_hPa'] = medfilt(df2['Pressure_hPa'], kernel_size=3)
                
                # Resample data: take the first point and then every 5th point
                df = df.iloc[::5].reset_index(drop=True)
                df2 = df2.iloc[::5].reset_index(drop=True)
                
                # Check for sufficient data
                if df.empty or df2.empty:
                    print(f"Skipping files {file} and {file2} due to insufficient data.")
                    continue
                
                # Convert to units
                p = df['Pressure_hPa'].values * units.hPa
                T = df['Temp_degC'].values * units.degC
                Td = df['Frostpoint_degC'].values * units.degC
                WindSpeed = df['WindSpeed_mps'].values * units.mps
                WindDir = df['WindDirection_deg'].values * units.deg
                
                p2 = df2['Pressure_hPa'].values * units.hPa
                T2 = df2['Temp_degC'].values * units.degC
                Td2 = df2['Frostpoint_degC'].values * units.degC
                WindSpeed2 = df2['WindSpeed_mps'].values * units.mps
                WindDir2 = df2['WindDirection_deg'].values * units.deg
                
                # Determine the u and v components of the wind
                u, v = mpcalc.wind_components(WindSpeed, WindDir)
                u2, v2 = mpcalc.wind_components(WindSpeed2, WindDir2)
                
                try:
                    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
                    lcl_pressure2, lcl_temperature2 = mpcalc.lcl(p2[0], T2[0], Td2[0])

                    el_pressure, el_temperature = mpcalc.el(p, T, Td)
                    el_pressure2, el_temperature2 = mpcalc.el(p2, T2, Td2)

                    lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td)
                    lfc_pressure2, lfc_temperature2 = mpcalc.lfc(p2, T2, Td2)

                    surface_cape, surface_cin = mpcalc.surface_based_cape_cin(p, T, Td)
                    surface_cape2, surface_cin2 = mpcalc.surface_based_cape_cin(p2, T2, Td2)

                    mu_cape, mu_cin = mpcalc.most_unstable_cape_cin(p, T, Td)
                    mu_cape2, mu_cin2 = mpcalc.most_unstable_cape_cin(p2, T2, Td2)

                    parcel_path = mpcalc.parcel_profile(p, T[0], Td[0])
                    parcel_path2 = mpcalc.parcel_profile(p2, T2[0], Td2[0])

                except Exception as e:
                    print(f"Skipping files {file} and {file2} due to error: {e}")
                    continue
                
                fig = plt.figure(figsize=(12,8))
        
                # SkewT
                    
                skew = plots.SkewT(fig)
                skew.ax.set_xlim(-90,30)
                skew.ax.set_ylim(1050,100)
                skew.plot(p, T, 'red', label="Temperature 1")
                skew.plot(p2, T2, 'red', linestyle="dashed", label="Temperature 2", alpha=0.6)

                skew.plot(p, Td, 'green', label="Dew Point 1")
                skew.plot(p2, Td2, 'green', linestyle="dashed", label="Dew Point 2", alpha=0.6)
                
                interval = np.logspace(2, 3) * units.hPa
                idx = mpcalc.resample_nn_1d(p, interval)
                idx2 = mpcalc.resample_nn_1d(p2, interval)

                skew.plot_barbs(p[idx][::2], u[idx][::2], v[idx][::2])
                skew.plot_barbs(p2[idx2][::2], u2[idx2][::2], v2[idx2][::2])

                
                skew.plot_dry_adiabats()
                skew.plot_moist_adiabats()
                skew.plot_mixing_lines()
                
                skew.plot(p, parcel_path, color='k')
                skew.plot(p2, parcel_path2, color='k', linestyle="dashed", alpha=0.6)

                skew.shade_cape(p, T, parcel_path)
                skew.shade_cin(p, T, parcel_path)
                
                if lcl_pressure:
                    skew.ax.axhline(lcl_pressure, color ='orange', label='LCL Pressure 1')
                if lcl_pressure2:
                    skew.ax.axhline(lcl_pressure2, color ='orange', label='LCL Pressure 2', linestyle="dashed", alpha=0.6)
                if lfc_pressure:
                    skew.ax.axhline(lfc_pressure, color='purple', label='LFC Pressure 1')
                if lfc_pressure2:
                    skew.ax.axhline(lfc_pressure2, color='purple', label='LFC Pressure 2', linestyle="dashed", alpha=0.6)
                if el_pressure:
                    skew.ax.axhline(el_pressure, color='blue', label='Equilibrium level 1')
                if el_pressure2:
                    skew.ax.axhline(el_pressure2, color='blue', label='Equilibrium level 2', linestyle="dashed", alpha=0.6) 
                
                skew.ax.legend(loc='upper right', bbox_to_anchor=(0.95, 1), fontsize=10, frameon=True)

            
                skew.ax.set_title('TRACER Sonde ' + launch_date_for_title, fontsize=10)
                fig_filename = os.path.join(output_folder, f"TRACERSonde.{launch_date}.png")
                      
                
                fig.savefig(fig_filename, bbox_inches='tight', dpi=300)
                plt.show()
            

combined_skewT(folder_path, output_folder)