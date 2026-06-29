import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np 
from scipy.optimize import curve_fit

#absorbance

ticks = [200, 300, 400, 500, 600, 700, 800, 900, 1000]

df_blue = pd.read_csv("./G5 ABS/G5/blue_absorption.csv", sep=";", header=32, names=["wavelength", "abs"])
df_blue = df_blue.dropna()
#df_blue["abs"] = -np.log10(df_blue["transmission"])
df_blue["wavelength"] = df_blue["wavelength"].astype(float)

df_green = pd.read_csv("./G5 ABS/G5/green_absorption.csv", sep=";", header=32, names=["wavelength", "abs"])
df_green = df_green.dropna()
#df_green["abs"] = -np.log10(df_green["transmission"])
df_green["wavelength"] = df_green["wavelength"].astype(float)

df_red = pd.read_csv("./G5 ABS/G5/red_absorption.csv", sep=";", header=32, names=["wavelength", "abs"])
df_red = df_red.dropna()
#df_red["abs"] = -np.log10(df_red["transmission"])
df_red["wavelength"] = df_red["wavelength"].astype(float)



plt.plot(df_blue["wavelength"], df_blue["abs"], label="Blue")
plt.plot(df_green["wavelength"], df_green["abs"], color="green", label="Green")
plt.plot(df_red["wavelength"], df_red["abs"], color ="red", label="Red")
plt.title("Absorbance spectrum ")
plt.legend()
plt.xticks(ticks)
plt.xlim((450,800))
plt.ylim((0, 4.5))
plt.xlabel("Wavelength [nm]")
plt.ylabel("Absorbance [a.u.]")
plt.savefig("absorbance_spec_overlaid.pdf")
plt.show()


#absorption coefficient

thickness_blue = 366.94*10**-7 # nm (10^-9) to cm (10^-2) factor 10^7
abs_coeff_blue = 2.303*df_blue["abs"]/thickness_blue

thickness_green = 331.04*10**-7 
abs_coeff_green = 2.303*df_green["abs"]/thickness_green

thickness_red = 680*10**-7
abs_coeff_red = 2.303*df_red["abs"]/thickness_red


df_red["wavelength"][1700:1900]
max(abs_coeff_red[1700:1900])
print(f"max_blue:{max(abs_coeff_blue[1000:1500])}, max green {max(abs_coeff_green[1180:1900])},max red {max(abs_coeff_red[1700:1900])} cm-1")

plt.plot(df_blue["wavelength"], abs_coeff_blue, label="Blue", color="blue") 
plt.title("Absorption coefficient")
plt.xticks(ticks)
plt.xlim((400,800))
plt.ylim((0, 100000))
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"$\alpha$ [cm$^{-1}$]")

plt.plot(df_green["wavelength"], abs_coeff_green, color="Green", label="Green") 
#plt.title("Absorption coefficient of green")
#plt.xticks(ticks)
#plt.xlim((400,800))
#plt.ylim((0, 100000))
#plt.xlabel("Wavelength [nm]")
#plt.ylabel(r"$\alpha$ [cm$^{-1}$]")
#plt.show()

plt.plot(df_red["wavelength"], abs_coeff_red, label="Red", color="red") 
plt.legend()
plt.savefig("Absorb_coeff_overlaid.pdf")
#plt.xticks(ticks)
#plt.xlim((400,800))
#plt.ylim((0, 100000))
#plt.xlabel("Wavelength [nm]")
#plt.ylabel(r"$\alpha$ [cm$^{-1}$]")
#plt.show()

#tauc plot

# first convert wavelength to eV
# this is equivalent to photon energy
def lin_func(x, a, b):
    return a*x + b

df_blue["eV"] = 1240/df_blue["wavelength"]
df_green["eV"] = 1240/df_green["wavelength"]
df_red["eV"] = 1240/df_red["wavelength"]

y_axis_tauc_blue = (abs_coeff_blue*df_blue["eV"])**2
x_axis_tauc_blue = df_blue["eV"]

y_axis_tauc_green = (abs_coeff_green*df_green["eV"])**2
x_axis_tauc_green = df_green["eV"]

y_axis_tauc_red = (abs_coeff_red*df_red["eV"])**2
x_axis_tauc_red = df_red["eV"]

# fit for the blue photodiode
popt_blue, pcov_blue = curve_fit(lin_func,x_axis_tauc_blue[1282:1324],y_axis_tauc_blue[1282:1324]/10**9)
x_plot_blue = np.linspace(2.5,2.6,1000)

plt.plot(x_axis_tauc_blue, y_axis_tauc_blue/10**9)
plt.title("Tauc plot blue", fontsize=20)
plt.xlim((2.4,2.8))
plt.ylim(0,50)
plt.xlabel(r"h$\nu$ [eV]", fontsize=18)
plt.ylabel(r"($\alpha$h$\nu$)$^2$/10$^9$ [(eV $\cdot$ cm$^{-1}$)$^2$]", fontsize=18)
#plt.plot(x_axis_tauc_blue[1282:1324],y_axis_tauc_blue[1282:1324]/10**9)
plt.plot(x_plot_blue, lin_func(x_plot_blue, *popt_blue))
plt.savefig("Tauc_plot_blue.pdf")
plt.show()

# fit for green
popt_green, pcov_green = curve_fit(lin_func,x_axis_tauc_green[1360:1447],y_axis_tauc_green[1360:1447]/10**9)
x_plot_green = np.linspace(2.2,2.37,1000)

plt.plot(x_axis_tauc_green, y_axis_tauc_green/10**9, color="green")
plt.title("Tauc plot green", fontsize=20)
plt.xlim((2.1,2.6))
plt.ylim(0,40)
plt.xlabel(r"h$\nu$ [eV]", fontsize=18)
plt.ylabel(r"($\alpha$h$\nu$)$^2$/10$^9$ [(eV $\cdot$ cm$^{-1}$)$^2$]", fontsize=18)
#plt.plot(x_axis_tauc_green[1360:1447], y_axis_tauc_green[1360:1447]/10**10)
plt.plot(x_plot_green, lin_func(x_plot_green, *popt_green), color="orange")
plt.savefig("Tauc_plot_green.pdf")
plt.show()
#fit for red 
popt_red, pcov_red = curve_fit(lin_func,x_axis_tauc_red[1690:1775],y_axis_tauc_red[1690:1775]/10**9)
x_plot_red = np.linspace(1.6,1.72,1000)


plt.plot(x_axis_tauc_red, y_axis_tauc_red/10**9, color="red")
plt.title("Tauc plot red", fontsize=20)
plt.xlim((1.5,1.75))
plt.xlabel(r"h$\nu$ [eV]", fontsize=18)
plt.ylabel(r"($\alpha$h$\nu$)$^2$/10$^9$ [(eV $\cdot$ cm$^{-1}$)$^2$]" , fontsize=18)
plt.ylim(0,10)
#plt.plot(x_axis_tauc_red[1690:1775], y_axis_tauc_red[1690:1775]/10**10)
plt.plot(x_plot_red, lin_func(x_plot_red, *popt_red), color="orange")
plt.savefig("Tauc_plot_red.pdf")
plt.show()

band_gap_blue = -popt_blue[1]/popt_blue[0]
nm_blue = 1240/band_gap_blue

band_gap_green = -popt_green[1]/popt_green[0]
nm_green = 1240/band_gap_green

band_gap_red = -popt_red[1]/popt_red[0]
nm_red = 1240/band_gap_red

print(f"The bandgaps are for blue:{band_gap_blue:.2f} eV, green:{band_gap_green:.2f} eV, red:{band_gap_red:.2f} eV")
print(f"blue: {nm_blue} green: {nm_green}, red: {nm_red} ")

