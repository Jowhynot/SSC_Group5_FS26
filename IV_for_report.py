import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
from scipy.optimize import curve_fit
import glob
import os 
from scipy.interpolate import interp1d

directory = "./"

file_pattern_dark = os.path.join(directory, "*DARK.dat")
all_files_dark = glob.glob(file_pattern_dark)

file_pattern_light = os.path.join(directory, "*LIGHT.dat")
all_files_light = glob.glob(file_pattern_light)
all_files_dark


light_dataframes = {}

for i in all_files_light:
    
    file_name = os.path.splitext(os.path.basename(i))[0]
    
    
    light_dataframes[file_name] = pd.read_csv(i, sep="\t", names=["V", "A"], header=0)
    
dark_dataframes = {} 

for i in all_files_dark: 
    file_name = os.path.splitext(os.path.basename(i))[0]
    
    dark_dataframes[file_name] = pd.read_csv(i, sep="\t", names=["V", "A"], header=0)
dark_dataframes

# IV curves
# Area of a pixel 
Area = 0.16 # cm^2 

# Blue 1
J_sc_blue = [] # short circuit current density --> current density at 0 V
V_oc_blue = [] # open circuit voltage --> voltage at which current 0 A
Responsivity_blue = []

for i in range(1, 7):
    
    key_light = f"BLUE1_{i}LIGHT" 
    key_dark = f"BLUE1_{i}DARK"
    
    if key_light in light_dataframes and key_dark in dark_dataframes:
        df_light = light_dataframes[key_light]
        df_dark = dark_dataframes[key_dark]
        title = f"1st Blue Pixel {i}"
        # IV curves 
        #plt.figure() 
        #plt.plot(df_light["V"], df_light["A"], label="illuminated")
        #plt.plot(df_dark["V"] ,df_dark["A"], label="dark")
        
        #plt.xlabel("Voltage [V]")
        #plt.ylabel("Current [A]")
        #plt.title(f"IV Curve for {title}")
        #plt.legend()
        #plt.savefig(directory+title+".pdf")
        #plt.show()
        
        # calculation of Responsivity at 0V
        P_opt = Area * 1 # optical power on the pixel in mW (LED intensity = 1 mW/cm^2)
        df_curr_light = df_light.sort_values(by="V").reset_index(drop=True)
        df_curr_dark = df_dark.sort_values(by="V").reset_index(drop=True)
        
        interp_light_func = interp1d(df_curr_light["V"], df_curr_light["A"]*1000, kind="linear")
        interp_dark_func = interp1d(df_curr_dark["V"], df_curr_dark["A"]*1000, kind="linear")
        
        I_light = interp_light_func(0.0)
        I_dark = interp_dark_func(0.0)
        
        
        I_photo = I_light - I_dark

        Responsivity_blue.append(I_photo/P_opt)
        
        
        # Jsc calculation 
        df_Jsc = df_light.sort_values(by="V").reset_index(drop=True) #sort dataframe by voltage
        interp_Jsc_func = interp1d(df_Jsc["V"], 1000*df_Jsc["A"]/Area, kind ="linear")
        
        val_Jsc = interp_Jsc_func(0.0)
        
        J_sc_blue.append(val_Jsc)
        
        # Voc calculation 
        df_Voc = df_light.sort_values(by="A").reset_index(drop=True) #sort dataframe by voltage
        interp_Voc_func = interp1d(df_Voc["A"]/Area, df_Voc["V"], kind ="linear")
        try:
            val_Voc = interp_Voc_func(0.0)
            V_oc_blue.append(val_Voc)
        except ValueError:
            print(r"No Voc found")
    
        # JV curves
        plt.figure()
        plt.plot(df_light["V"], 1000*df_light["A"]/Area, label="illuminated")
        plt.plot(df_dark["V"] ,1000*df_dark["A"]/Area, label="dark")
        plt.xlabel("Voltage [V]")
        plt.ylabel(r"Current density [mA/cm$^2$]")
        plt.title(f"JV Curve for {title}")
        #plt.xlim((-0.01, 0.01))
        plt.legend()
        plt.savefig(directory+"JV"+title+".pdf")
        plt.show()
    else:
        print(f"Warning: {key} not found in dataframes.")
    
# Blue 2

for i in range(1, 7):
    
    key_light = f"BLUE2_{i}LIGHT" 
    key_dark = f"BLUE2_{i}DARK"
    
    if key_light in light_dataframes and key_dark in dark_dataframes:
        df_light = light_dataframes[key_light]
        df_dark = dark_dataframes[key_dark]
        title = f"2nd Blue Pixel {i}"
        
        #plt.figure() 
        #plt.plot(df_light["V"], df_light["A"], label="illuminated")
        #plt.plot(df_dark["V"] ,df_dark["A"], label="dark")
        #
        #plt.xlabel("Voltage [V]")
        #plt.ylabel("Current [A]")
        #plt.title(f"IV Curve for {title}")
        #plt.legend()
        #plt.savefig(directory+title+".pdf")
        #plt.show()
        
       # calculation of Responsivity at 0V
        P_opt = Area * 1 # optical power on the pixel in mW (LED intensity = 1 mW/cm^2)
        df_curr_light = df_light.sort_values(by="V").reset_index(drop=True)
        df_curr_dark = df_dark.sort_values(by="V").reset_index(drop=True)
        
        interp_light_func = interp1d(df_curr_light["V"], df_curr_light["A"]*1000, kind="linear")
        interp_dark_func = interp1d(df_curr_dark["V"], df_curr_dark["A"]*1000, kind="linear")
        
        I_light = interp_light_func(0.0)
        I_dark = interp_dark_func(0.0)
        
        
        I_photo = I_light - I_dark

        Responsivity_blue.append(I_photo/P_opt)
        
        #Jsc calculation 
        df_Jsc = df_light.sort_values(by="V").reset_index(drop=True) #sort dataframe by voltage
        interp_Jsc_func = interp1d(df_Jsc["V"], 1000*df_Jsc["A"]/Area, kind ="linear")
        
        val_Jsc = interp_Jsc_func(0.0)
        
        J_sc_blue.append(val_Jsc)
        
        # Voc calculation 
        df_Voc = df_light.sort_values(by="A").reset_index(drop=True) #sort dataframe by voltage
        interp_Voc_func = interp1d(df_Voc["A"]/Area, df_Voc["V"], kind ="linear")
        
        try: 
            val_Voc = interp_Voc_func(0.0)
            V_oc_blue.append(val_Voc)
            
        except ValueError: 
            print("No Voc found")
        
        # JV curves
        plt.figure()
        plt.plot(df_light["V"], 1000*df_light["A"]/Area, label="illuminated")
        plt.plot(df_dark["V"] ,1000*df_dark["A"]/Area, label="dark")
        plt.xlabel("Voltage [V]")
        plt.ylabel(r"Current density [mA/cm$^2$]")
        plt.title(f"JV Curve for {title}")
        #plt.xlim((-0.001, 0.001))
        #plt.ylim((-3,0.1))
        plt.legend()
        plt.savefig(directory+"JV"+title+".pdf")
        plt.show()
    else:
        print(f"Warning: {key} not found in dataframes.")

Responsivity_blue = np.array(Responsivity_blue)

# Green 1
J_sc_green = [] # short circuit current density --> current density at 0 V
V_oc_green = [] # open circuit voltage --> voltage at which current 0 A
Responsivity_green = []

for i in range(1, 7):
    
    key_light = f"GREEN1_{i}LIGHT" 
    key_dark = f"GREEN1_{i}DARK"
    
    if key_light in light_dataframes and key_dark in dark_dataframes:
        df_light = light_dataframes[key_light]
        df_dark = dark_dataframes[key_dark]
        title = f"1st Green Pixel {i}"
        
        #plt.figure() 
        #plt.plot(df_light["V"], df_light["A"], label="illuminated")
        #plt.plot(df_dark["V"] ,df_dark["A"], label="dark")
        #
        #plt.xlabel("Voltage [V]")
        #plt.ylabel("Current [A]")
        #plt.title(f"IV Curve for {title}")
        #plt.legend()
        #plt.savefig(directory+title+".pdf")
        #plt.show()
        
        # calculation of Responsivity at 0V
        P_opt = Area * 1 # optical power on the pixel in mW (LED intensity = 1 mW/cm^2)
        df_curr_light = df_light.sort_values(by="V").reset_index(drop=True)
        df_curr_dark = df_dark.sort_values(by="V").reset_index(drop=True)
        
        interp_light_func = interp1d(df_curr_light["V"], df_curr_light["A"]*1000, kind="linear")
        interp_dark_func = interp1d(df_curr_dark["V"], df_curr_dark["A"]*1000, kind="linear")
        
        I_light = interp_light_func(0.0)
        I_dark = interp_dark_func(0.0)
        
        
        I_photo = I_light - I_dark

        Responsivity_green.append(I_photo/P_opt)
        
        #J_sc calculation 
        df_Jsc = df_light.sort_values(by="V").reset_index(drop=True) #sort dataframe by voltage
        interp_Jsc_func = interp1d(df_Jsc["V"], 1000*df_Jsc["A"]/Area, kind ="linear")
        
        val_Jsc = interp_Jsc_func(0.0)
        
        J_sc_green.append(val_Jsc)
        
        # Voc calculation 
        df_Voc = df_light.sort_values(by="A").reset_index(drop=True) #sort dataframe by voltage
        interp_Voc_func = interp1d(df_Voc["A"]/Area, df_Voc["V"], kind ="linear")
        
        val_Voc = interp_Voc_func(0.0)
        
        V_oc_green.append(val_Voc)

        # JV curves
        plt.figure()
        plt.plot(df_light["V"], 1000*df_light["A"]/Area, label="illuminated")
        plt.plot(df_dark["V"] ,1000*df_dark["A"]/Area, label="dark")
        plt.xlabel("Voltage [V]")
        plt.ylabel(r"Current density [mA/cm$^2$]")
        plt.title(f"JV Curve for {title}")
        plt.legend()
        plt.savefig(directory+"JV"+title+".pdf")
        plt.show()
    else:
        print(f"Warning: {key} not found in dataframes.")
        
    # Green 2
for i in range(1, 7):
    
    key_light = f"GREEN2_{i}LIGHT" 
    key_dark = f"GREEN2_{i}DARK"
    
    if key_light in light_dataframes and key_dark in dark_dataframes:
        df_light = light_dataframes[key_light]
        df_dark = dark_dataframes[key_dark]
        title = f"2nd Green Pixel {i}"
        
        
        # calculation of Responsivity at 0V
        P_opt = Area * 1 # optical power on the pixel in mW (LED intensity = 1 mW/cm^2)
        df_curr_light = df_light.sort_values(by="V").reset_index(drop=True)
        df_curr_dark = df_dark.sort_values(by="V").reset_index(drop=True)
        
        interp_light_func = interp1d(df_curr_light["V"], df_curr_light["A"]*1000, kind="linear")
        interp_dark_func = interp1d(df_curr_dark["V"], df_curr_dark["A"]*1000, kind="linear")
        
        I_light = interp_light_func(0.0)
        I_dark = interp_dark_func(0.0)
        
        
        I_photo = I_light - I_dark

        Responsivity_green.append(I_photo/P_opt)
        
        #J_sc calculation 
        df_Jsc = df_light.sort_values(by="V").reset_index(drop=True) #sort dataframe by voltage
        interp_Jsc_func = interp1d(df_Jsc["V"], 1000*df_Jsc["A"]/Area, kind ="linear")
        
        val_Jsc = interp_Jsc_func(0.0)
        
        J_sc_green.append(val_Jsc)
        
        # Voc calculation 
        df_Voc = df_light.sort_values(by="A").reset_index(drop=True) #sort dataframe by voltage
        interp_Voc_func = interp1d(df_Voc["A"]/Area, df_Voc["V"], kind ="linear")
        
        val_Voc = interp_Voc_func(0.0)
        
        V_oc_green.append(val_Voc)
        
        # JV curves
        plt.figure()
        plt.plot(df_light["V"], 1000*df_light["A"]/Area, label="illuminated")
        plt.plot(df_dark["V"] ,1000*df_dark["A"]/Area, label="dark")
        plt.xlabel("Voltage [V]")
        plt.ylabel(r"Current density [mA/cm$^2$]")
        plt.title(f"JV Curve for {title}")
        plt.legend()
        plt.savefig(directory+"JV"+title+".pdf")
        plt.show()
    else:
        print(f"Warning: {key} not found in dataframes.")
        
Responsivity_green = np.array(Responsivity_green)

# Red 1
J_sc_red = []
V_oc_red = [] 
Responsivity_red = []

for i in range(1, 7):
    
    key_light = f"RED1_{i}LIGHT" 
    key_dark = f"RED1_{i}DARK"
    
    if key_light in light_dataframes and key_dark in dark_dataframes:
        df_light = light_dataframes[key_light]
        df_dark = dark_dataframes[key_dark]
        title = f"1st Red Pixel {i}"
        
        #plt.figure() 
        #plt.plot(df_light["V"], df_light["A"], label="illuminated")
        #plt.plot(df_dark["V"] ,df_dark["A"], label="dark")
        #
        #plt.xlabel("Voltage [V]")
        #plt.ylabel("Current [A]")
        #plt.title(f"IV Curve for {title}")
        #plt.legend()
        #plt.savefig(directory+title+".pdf")
        #plt.show()
        
        # calculation of Responsivity at 0V
        P_opt = Area * 1 # optical power on the pixel in mW (LED intensity = 1 mW/cm^2)
        df_curr_light = df_light.sort_values(by="V").reset_index(drop=True)
        df_curr_dark = df_dark.sort_values(by="V").reset_index(drop=True)
        
        interp_light_func = interp1d(df_curr_light["V"], df_curr_light["A"]*1000, kind="linear")
        interp_dark_func = interp1d(df_curr_dark["V"], df_curr_dark["A"]*1000, kind="linear")
        
        I_light = interp_light_func(0.0)
        I_dark = interp_dark_func(0.0)
        
        
        I_photo = I_light - I_dark

        Responsivity_red.append(I_photo/P_opt)
        
        # Jsc calculation 
        df_Jsc = df_light.sort_values(by="V").reset_index(drop=True) #sort dataframe by voltage
        interp_Jsc_func = interp1d(df_Jsc["V"], 1000*df_Jsc["A"]/Area, kind ="linear")
        
        val_Jsc = interp_Jsc_func(0.0)
        
        J_sc_red.append(val_Jsc)
        
        # Voc calculation 
        df_Voc = df_light.sort_values(by="A").reset_index(drop=True) #sort dataframe by voltage
        interp_Voc_func = interp1d(df_Voc["A"]/Area, df_Voc["V"], kind ="linear")
        
        val_Voc = interp_Voc_func(0.0)
        
        V_oc_red.append(val_Voc)
        
        # JV curves
        plt.figure()
        plt.plot(df_light["V"], 1000*df_light["A"]/Area, label="illuminated")
        plt.plot(df_dark["V"] ,1000*df_dark["A"]/Area, label="dark")
        plt.xlabel("Voltage [V]")
        plt.ylabel(r"Current density [mA/cm$^2$]")
        plt.title(f"JV Curve for {title}")
        plt.legend()
        plt.savefig(directory+"JV"+title+".pdf")
        plt.show()
    else:
        print(f"Warning: {key} not found in dataframes.")
        
# Red 2

for i in range(1, 7):
    
    key_light = f"RED2_{i}LIGHT" 
    key_dark = f"RED2_{i}DARK"
    
    if key_light in light_dataframes and key_dark in dark_dataframes:
        df_light = light_dataframes[key_light]
        df_dark = dark_dataframes[key_dark]
        title = f"2nd Red Pixel {i}"
        
        #plt.figure() 
        #plt.plot(df_light["V"], df_light["A"], label="illuminated")
        #plt.plot(df_dark["V"] ,df_dark["A"], label="dark")
        #
        #plt.xlabel("Voltage [V]")
        #plt.ylabel("Current [A]")
        #plt.title(f"IV Curve for {title}")
        #plt.legend()
        #plt.savefig(directory+title+".pdf")
        #plt.show()
        
        # calculation of Responsivity at 0V
        P_opt = Area * 1 # optical power on the pixel in mW (LED intensity = 1 mW/cm^2)
        df_curr_light = df_light.sort_values(by="V").reset_index(drop=True)
        df_curr_dark = df_dark.sort_values(by="V").reset_index(drop=True)
        
        interp_light_func = interp1d(df_curr_light["V"], df_curr_light["A"]*1000, kind="linear")
        interp_dark_func = interp1d(df_curr_dark["V"], df_curr_dark["A"]*1000, kind="linear")
        
        I_light = interp_light_func(0.0)
        I_dark = interp_dark_func(0.0)
        
        
        I_photo = I_light - I_dark

        Responsivity_red.append(I_photo/P_opt)
        
        # Jsc calculation 
        df_Jsc = df_light.sort_values(by="V").reset_index(drop=True) #sort dataframe by voltage
        interp_Jsc_func = interp1d(df_Jsc["V"], 1000*df_Jsc["A"]/Area, kind ="linear")
        
        val_Jsc = interp_Jsc_func(0.0)
        
        J_sc_red.append(val_Jsc)
        
        # Voc calculation 
        df_Voc = df_light.sort_values(by="A").reset_index(drop=True) #sort dataframe by voltage
        interp_Voc_func = interp1d(df_Voc["A"]/Area, df_Voc["V"], kind ="linear")
        
        val_Voc = interp_Voc_func(0.0)
        
        V_oc_red.append(val_Voc)
        
        # JV curves
        plt.figure()
        plt.plot(df_light["V"], 1000*df_light["A"]/Area, label="illuminated")
        plt.plot(df_dark["V"] ,1000*df_dark["A"]/Area, label="dark")
        plt.xlabel("Voltage [V]")
        plt.ylabel(r"Current density [mA/cm$^2$]")
        plt.title(f"JV Curve for {title}")
        plt.legend()
        plt.savefig(directory+"JV"+title+".pdf")
        plt.show()
        
    else:
        print(f"Warning: {key} not found in dataframes.")
        
Responsivity_red = np.array(Responsivity_red)

# results 

# V_oc and Js

print(f"The Jsc of the blue pixels: {np.mean(J_sc_blue):.4f} \u00B1 {np.std(J_sc_blue):.4f}mA/cm^2, green pixels: 
      {np.mean(J_sc_green):.4f} \u00B1 {np.std(J_sc_green):.4f}mA/cm^2, red pixels: 
      {np.mean(J_sc_red):.4f} \u00B1 {np.std(J_sc_red):.4f}mA/cm^2")

print(f"The Voc of the blue pixels: {np.mean(V_oc_blue):.2f} \u00B1 {np.std(V_oc_blue):.2f}V, green pixels: 
      {np.mean(V_oc_green):.2f} \u00B1 {np.std(V_oc_green):.2f}V, 
      red pixels: {np.mean(V_oc_red):.2f} \u00B1 {np.std(V_oc_red) :.2f}V")

# Responsivity 
print(f"The responsivity of the blue pixels: {np.mean(Responsivity_blue):.2f} \u00B1 {np.std(Responsivity_blue):.2f}A/W, 
      green pixels: {np.mean(Responsivity_green):.2f} \u00B1 {np.std(Responsivity_green):.2f}A/W, 
      red pixels: {np.mean(Responsivity_red):.2f} \u00B1 {np.std(Responsivity_red):.2f}A/W")  
 
# EQE 
wavelength_LED_blue = 461*10**-9# 
wavelength_LED_green = 522*10**-9 #
wavelength_LED_red = 626*10**-9 #
h = 6.6261 * 10**-34# plancks constant J/s
c = 299792458 # light speed m/s
q = 1.602*10**-19 # elementary charge C 
EQE_blue = (Responsivity_blue*h*c)/(q*wavelength_LED_blue)*100
EQE_green = (Responsivity_green*h*c)/(q*wavelength_LED_green)*100
EQE_red = (Responsivity_red*h*c)/(q*wavelength_LED_red)*100

print(f"The EQE of blue pixels is: {np.mean(EQE_blue):.2f} \u00B1 {np.std(EQE_blue):.2f}%, 
      green pixels: {np.mean(EQE_green):.2f} \u00B1 {np.std(EQE_green):.2f}%, 
      red pixels: {np.mean(EQE_red):.2f} \u00B1 {np.std(EQE_red):.2f}%")


#J_sc
plt.bar("Blue",np.mean(J_sc_blue), yerr=np.std(J_sc_blue), color ="Blue")
plt.bar("Green", np.mean(J_sc_green), yerr=np.std(J_sc_green), color="Green")
plt.bar("Red", np.mean(J_sc_red), yerr=np.std(J_sc_red), color="Red")
plt.ylabel(r"J$_{sc}$ [mA/cm$^2$]")
plt.savefig("err_plot_Jsc.pdf")
plt.show()

#V_oc 
plt.bar("Blue",np.mean(V_oc_blue), yerr=np.std(V_oc_blue), color ="Blue")
plt.bar("Green", np.mean(V_oc_green), yerr=np.std(V_oc_green), color="Green")
plt.bar("Red", np.mean(V_oc_red), yerr=np.std(V_oc_red), color="Red")
plt.ylabel(r"V$_{oc}$ [V]")
plt.savefig("err_plot_Voc.pdf")
plt.show()

# Responsivity 
plt.bar("Blue",np.abs(np.mean(Responsivity_blue)), yerr=np.std(Responsivity_blue), color ="Blue")
plt.bar("Green", np.abs(np.mean(Responsivity_green)), yerr=np.std(Responsivity_green), color="Green")
plt.bar("Red", np.abs(np.mean(Responsivity_red)), yerr=np.std(Responsivity_red), color="Red")
plt.ylabel(r"Responsivity [A/W]")
plt.savefig("err_plot_responsivity.pdf")
plt.show()

# EQE 
plt.bar("Blue",np.abs(np.mean(EQE_blue)), yerr=np.abs(np.std(EQE_blue)), color ="Blue")
plt.bar("Green", np.abs(np.mean(EQE_green)), yerr=np.abs(np.std(EQE_green)), color="Green")
plt.bar("Red", np.abs(np.mean(EQE_red)), yerr=np.abs(np.std(EQE_red)), color="Red")
plt.ylabel(r"EQE [%]")
plt.savefig("err_plot_eqe.pdf")
plt.show()