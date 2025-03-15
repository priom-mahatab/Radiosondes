#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:59:22 2025

@author: zmahatab
"""
import os
from datetime import datetime
from PIL import Image, ImageChops

def trim_whitespace(image_path, padding=300):
    img = Image.open(image_path).convert("RGB")
    bg = Image.new("RGB", img.size, (255,255,255))
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    
    if bbox:
        x1, y1, x2, y2 = bbox
        x1 = max(0, x1 - padding)
        y1 = max(0, y1 - padding)
        x2 = min(img.width, x2 + padding)
        y2 = min(img.height, y2 + padding)
        return img.crop((x1, y1, x2, y2))
        
    return img

input_folder = "/Users/zmahatab/Desktop/MetPy/TRACER_M1/Plots"
input_folder_2 = "/Users/zmahatab/Desktop/MetPy/TRACER_M1/Plots/Combined"
output_folder = "/Users/zmahatab/Desktop/MetPy/TRACER_M1/Plots/Three-Grid"


def three_grid_generator(input_folder, input_folder_2, output_folder):
    for file1 in os.listdir(input_folder):
        if file1.endswith(".png"):
            parts1 = file1.split(".")
            if len(parts1) == 5:
                date_part1 = parts1[1] # YYYYMMDD
                time_part1 = parts1[2] # HHMMSS
                hour_part1 = parts1[2][:2] # HH
                
                if hour_part1 == "11":
                    file_path1 = os.path.join(input_folder, file1)
                    
                    formatted_date = f"{date_part1[4:6]}-{date_part1[6:8]}-{date_part1[:4]}"
                    
                    print(f"Processing {file1}")
                    
                    trimmed_image1 = trim_whitespace(file_path1)
                    
                    for file2 in os.listdir(input_folder):
                        if file2.endswith(".png"):
                            parts2 = file2.split(".")
                            date_part2 = parts2[1]
                            time_part2 = parts2[2][:2] ## HH
                            
                            if date_part1 == date_part2 and time_part2 == "17":
                                file_path2 = os.path.join(input_folder, file2)
                                print(f"Processing {file2}")
                                
                                trimmed_image2 = trim_whitespace(file_path2)
                    
                    for file3 in os.listdir(input_folder_2):
                        if file3.endswith(".png"):
                            parts3 = file3.split(".")
                            date_part3 = parts3[1]
                            if date_part1 == date_part3:
                                file_path3 = os.path.join(input_folder_2, file3)
                                print(f"Processing {file3}")
                                trimmed_image3 = trim_whitespace(file_path3)
                                
            
                                height = min(trimmed_image1.height, trimmed_image2.height, trimmed_image3.height)
                                width = trimmed_image1.width + trimmed_image2.width + trimmed_image3.width
                                combined_image = Image.new('RGB', (width, height))
                                    
                                combined_image.paste(trimmed_image1, (0,0))
                                combined_image.paste(trimmed_image2, (trimmed_image1.width, 0))
                                combined_image.paste(trimmed_image3, (trimmed_image1.width + trimmed_image2.width, 0))                                
                                    
                                #combined_image.show()
                                
                                output_filename = os.path.join(output_folder, f"TRACERM1.{formatted_date}.three-grid.png")
                                
                                combined_image.save(output_filename)


three_grid_generator(input_folder, input_folder_2, output_folder)
                    