import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

#IV curves 
#gold gold contacted device 
df_neg1 = pd.read_csv("Au_Au_neg1V_79uA_V2.csv")
df_neg3 = pd.read_csv("Au_Au_neg3V_79uA.csv")
df_neg5 = pd.read_csv("Au_Au_neg5V_79uA_V2.csv")
df_neg7 = pd.read_csv("Au_Au_neg7V_79uA.csv")
df_neg9 = pd.read_csv("Au_Au_neg9V_79uA.csv")

# get the current from the equation U/R=I 
R = 10**6 # 1 Mega Ohm

df_neg1["Current A"] = df_neg1["Vsense 1 V (V)"]/R
df_neg3["Current A"] = df_neg3["Vsense 1 V (V)"]/R
df_neg5["Current A"] = df_neg5["Vsense 1 V (V)"]/R
df_neg7["Current A"] = df_neg7["Vsense 1 V (V)"]/R
df_neg9["Current A"] = df_neg9["Vsense 1 V (V)"]/R


plt.plot(df_neg1["Current A"])
plt.title(r"Gold Gold -1V 80 $\mu$A")
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
#plt.plot(df_neg1["Current A"].iloc[750:820])
plt.savefig("AuAu_xrays_neg1V.pdf")
plt.savefig("AuAu_xrays_neg1V.eps")
plt.show()



plt.plot(df_neg3["Current A"])
plt.title(r"Gold Gold -3V 80 $\mu$A")
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
#plt.plot(df_neg3["Current A"].iloc[1800:1980])
plt.savefig("AuAu_xrays_neg3V.pdf")
plt.savefig("AuAu_xrays_neg3V.eps")
plt.show()



plt.plot(df_neg5["Current A"])
plt.title(r"Gold Gold -5V 80 $\mu$A")
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
#plt.plot(df_neg5["Current A"].iloc[675:765])
plt.savefig("AuAu_xrays_neg5V.pdf")
plt.savefig("AuAu_xrays_neg5V.eps")
plt.show()




plt.plot(df_neg7["Current A"])
plt.title(r"Gold Gold -7V 80 $\mu$A")
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
#plt.plot(df_neg7["Current A"].iloc[895:980])
plt.savefig("AuAu_xrays_neg7V.pdf")
plt.savefig("AuAu_xrays_neg7V.eps")
plt.show()



plt.plot(df_neg9["Current A"])
plt.title(r"Gold Gold -9V 80 $\mu$A")
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
#plt.plot(df_neg9["Current A"].iloc[405:468])
plt.savefig("AuAu_xrays_neg9V.pdf")
plt.savefig("AuAu_xrays_neg9V.eps")
plt.show()

curr_neg9 = []

neg9_dark1 = np.mean(df_neg9["Current A"].iloc[20:70])
neg9_light1 = np.mean(df_neg9["Current A"].iloc[80:130])
curr_neg9.append(neg9_light1-neg9_dark1) 

neg9_dark2 = np.mean(df_neg9["Current A"].iloc[140:220])
neg9_light2 = np.mean(df_neg9["Current A"].iloc[230:295])
curr_neg9.append(neg9_light2-neg9_dark2) 

neg9_dark3 = np.mean(df_neg9["Current A"].iloc[305:395])
neg9_light3 = np.mean(df_neg9["Current A"].iloc[405:468])
curr_neg9.append(neg9_light3-neg9_dark3) 


#Hecht model definition
# Hecht Model definiton 
d_crystal = 0.00107 # thickness of the crystal m 
d_molfilm = 0.000042 # thickness of the molten film m 
def Hecht_eq_crystal(V, mutau, I_0): 
    d=d_crystal
    return I_0*mutau*V/(d**2)*(1-np.exp(-(d**2)/(mutau*V)))
def Hecht_eq_molfilm(V, mutau, I_0): 
    d=d_molfilm
    return I_0*mutau*V/(d**2)*(1-np.exp(-(d**2)/(mutau*V)))


curr_neg1 = []

neg1_dark1 = np.mean(df_neg1["Current A"].iloc[0:50])
neg1_light1 = np.mean(df_neg1["Current A"].iloc[70:160])
curr_neg1.append(neg1_light1-neg1_dark1)

neg1_dark2 = np.mean(df_neg1["Current A"].iloc[180:325])
neg1_light2 = np.mean(df_neg1["Current A"].iloc[335:430])
curr_neg1.append(neg1_light2-neg1_dark2)

neg1_dark3 = np.mean(df_neg1["Current A"].iloc[450:555])
neg1_light3 = np.mean(df_neg1["Current A"].iloc[565:630])
curr_neg1.append(neg1_light3-neg1_dark3)

neg1_dark4 = np.mean(df_neg1["Current A"].iloc[645:745])
neg1_light4 = np.mean(df_neg1["Current A"].iloc[750:820])
curr_neg1.append(neg1_light4-neg1_dark4)


 

curr_neg3 = []

neg3_dark1 = np.mean(df_neg3["Current A"].iloc[10:310])
neg3_light1 = np.mean(df_neg3["Current A"].iloc[370:570])
curr_neg3.append(neg3_light1-neg3_dark1)

neg3_dark2 = np.mean(df_neg3["Current A"].iloc[650:950])
neg3_light2 = np.mean(df_neg3["Current A"].iloc[1000:1280])
curr_neg3.append(neg3_light2-neg3_dark2) 

neg3_dark3 = np.mean(df_neg3["Current A"].iloc[1330:1730])
neg3_light3 = np.mean(df_neg3["Current A"].iloc[1800:1980])
curr_neg3.append(neg3_light3-neg3_dark3) 



curr_neg5 = []

neg5_dark1 = np.mean(df_neg5["Current A"].iloc[0:50])
neg5_light1 = np.mean(df_neg5["Current A"].iloc[80:150])
curr_neg5.append(neg5_light1-neg5_dark1)

neg5_dark2 = np.mean(df_neg5["Current A"].iloc[160:260])
neg5_light2 = np.mean(df_neg5["Current A"].iloc[270:340])
curr_neg5.append(neg5_light2-neg5_dark2)

neg5_dark3 = np.mean(df_neg5["Current A"].iloc[350:450])
neg5_light3 = np.mean(df_neg5["Current A"].iloc[465:530])
curr_neg5.append(neg5_light3-neg5_dark3)

neg5_dark4 = np.mean(df_neg5["Current A"].iloc[545:670])
neg5_light4 = np.mean(df_neg5["Current A"].iloc[675:765])
curr_neg5.append(neg5_light4-neg5_dark4)



curr_neg7 = []

neg7_dark1 = np.mean(df_neg7["Current A"].iloc[0:100])
neg7_light1 = np.mean(df_neg7["Current A"].iloc[110:155])
curr_neg7.append(neg7_light1-neg7_dark1) 

neg7_dark2 = np.mean(df_neg7["Current A"].iloc[180:330])
neg7_light2 = np.mean(df_neg7["Current A"].iloc[340:431])
curr_neg7.append(neg7_light2-neg7_dark2) 

neg7_dark3 = np.mean(df_neg7["Current A"].iloc[445:570])
neg7_light3 = np.mean(df_neg7["Current A"].iloc[580:639])
curr_neg7.append(neg7_light3-neg7_dark3) 

neg7_dark4 = np.mean(df_neg7["Current A"].iloc[650:733])
neg7_light4 = np.mean(df_neg7["Current A"].iloc[740:820])
curr_neg7.append(neg7_light4-neg7_dark4) 

neg7_dark5 = np.mean(df_neg7["Current A"].iloc[830:883])
neg7_light5 = np.mean(df_neg7["Current A"].iloc[895:980])
curr_neg7.append(neg7_light5-neg7_dark5)




curr_neg9 = []

neg9_dark1 = np.mean(df_neg9["Current A"].iloc[20:70])
neg9_light1 = np.mean(df_neg9["Current A"].iloc[80:130])
curr_neg9.append(neg9_light1-neg9_dark1) 

neg9_dark2 = np.mean(df_neg9["Current A"].iloc[140:220])
neg9_light2 = np.mean(df_neg9["Current A"].iloc[230:295])
curr_neg9.append(neg9_light2-neg9_dark2) 

neg9_dark3 = np.mean(df_neg9["Current A"].iloc[305:395])
neg9_light3 = np.mean(df_neg9["Current A"].iloc[405:468])
curr_neg9.append(neg9_light3-neg9_dark3) 

photo_current_au = np.abs([np.mean(curr_neg1), np.mean(curr_neg3), np.mean(curr_neg5), np.mean(curr_neg7), np.mean(curr_neg9)])

bias_gold_gold = np.abs([-1, -3, -5, -7, -9])
plt.scatter(bias_gold_gold, photo_current_au*10**9)
plt.ylabel("Photocurrent [pA]")
plt.xlabel("Absolute Voltage [V]")

# fit 
popt_au_au, pcov_au_au = curve_fit(Hecht_eq_crystal, bias_gold_gold, photo_current_au, p0=[0.000001, 200])
x_plot_au = np.linspace(0.05,10, 1000)
plt.plot(x_plot_au, Hecht_eq_crystal(x_plot_au, *popt_au_au)*10**9, color="coral")
plt.savefig("Hecht_fit_auau.pdf")
print(f"mobility lifetime product is {popt_au_au[0]}m^2/V with standard error {np.sqrt(pcov_au_au[0][0])}m^2/V")

#Gold Bismuth 

#df_bi_10V = pd.read_csv("Au_Bi_10V_79uA.csv") #not used
#df_bi_15V = pd.read_csv("Au_Bi_15V_79uA.csv") # not used 
df_bi_21V = pd.read_csv("Au_Bi_21V_80uA_V3.csv")
df_bi_41V = pd.read_csv("Au_Bi_41V_80uA_V3.csv")
df_bi_81V = pd.read_csv("Au_Bi_81V_80uA_V2.csv")
df_bi_101V = pd.read_csv("Au_Bi_101V_80uA.csv")


#df_bi_10V["Current A"] = df_bi_10V["Vsense 1 V (V)"]/R
#df_bi_15V["Current A"] = df_bi_15V["Vsense 1 V (V)"]/R
df_bi_21V["Current A"] = df_bi_21V["Vsense 1 V (V)"]/R
df_bi_41V["Current A"] = df_bi_41V["Vsense 1 V (V)"]/R
df_bi_81V["Current A"] = df_bi_81V["Vsense 1 V (V)"]/R
df_bi_101V["Current A"] = df_bi_101V["Vsense 1 V (V)"]/R

plt.plot(df_bi_21V["Current A"])
plt.title(r"Gold Bismuth 21V 80 $\mu$A")
#plt.plot(df_bi_21V["Current A"].iloc[500:596])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("AuBi_xrays_21V.pdf")
plt.savefig("AuBi_xrays_21V.eps")
plt.show()

plt.plot(df_bi_41V["Current A"])
plt.title(r"Gold Bismuth 41V 80 $\mu$A")
#plt.plot(df_bi_41V["Current A"].iloc[593:710])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("AuBi_xrays_41V.pdf")
plt.savefig("AuBi_xrays_41V.eps")
plt.show()

plt.plot(df_bi_81V["Current A"])
plt.title(r"Gold Bismuth 81V 80 $\mu$A")
#plt.plot(df_bi_81V["Current A"].iloc[960:1100])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("AuBi_xrays_81V.pdf")
plt.savefig("AuBi_xrays_81V.eps")
plt.show()

plt.plot(df_bi_101V["Current A"])
plt.title(r"Gold Bismuth 101V 80 $\mu$A")
#plt.plot(df_bi_101V["Current A"].iloc[1100:1196])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("AuBi_xrays_101V.pdf")
plt.savefig("AuBi_xrays_101V.eps")
plt.show()

curr_bi_21 = []

V21_dark1 = np.mean(df_bi_21V["Current A"].iloc[0:44])
V21_light1 = np.mean(df_bi_21V["Current A"].iloc[46:145])
curr_bi_21.append(V21_light1-V21_dark1) 

V21_dark2 = np.mean(df_bi_21V["Current A"].iloc[150:270])
V21_light2 = np.mean(df_bi_21V["Current A"].iloc[275:355])
curr_bi_21.append(V21_light2-V21_dark2)

V21_dark3 = np.mean(df_bi_21V["Current A"].iloc[360:490])
V21_light3 = np.mean(df_bi_21V["Current A"].iloc[500:596])
curr_bi_21.append(V21_light3-V21_dark3)

curr_bi_41 = []

V41_dark1 = np.mean(df_bi_41V["Current A"].iloc[0:135])
V41_light1 = np.mean(df_bi_41V["Current A"].iloc[146:305])
curr_bi_41.append(V41_light1-V41_dark1)

V41_dark2 = np.mean(df_bi_41V["Current A"].iloc[310:433])
V41_light2 = np.mean(df_bi_41V["Current A"].iloc[440:587])
curr_bi_41.append(V41_light2-V41_dark2)

V41_dark3 = np.mean(df_bi_41V["Current A"].iloc[593:710])
V41_light3 = np.mean(df_bi_41V["Current A"].iloc[720:833])
curr_bi_41.append(V41_light3-V41_dark3)

curr_bi_81 = []

V81_dark1 = np.mean(df_bi_81V["Current A"].iloc[150:251])
V81_light1 = np.mean(df_bi_81V["Current A"].iloc[260:430])
curr_bi_81.append(V81_light1-V81_dark1) 

V81_dark2 = np.mean(df_bi_81V["Current A"].iloc[470:630])
V81_light2 = np.mean(df_bi_81V["Current A"].iloc[640:800])
curr_bi_81.append(V81_light2-V81_dark2) 

V81_dark3 = np.mean(df_bi_81V["Current A"].iloc[840:940])
V81_light3 = np.mean(df_bi_81V["Current A"].iloc[960:1100])
curr_bi_81.append(V81_light3-V81_dark3) 

curr_bi_101 = []

V101_dark1 = np.mean(df_bi_101V["Current A"].iloc[0:80])
V101_light1 = np.mean(df_bi_101V["Current A"].iloc[90:140])
curr_bi_101.append(V101_light1-V101_dark1) 

V101_dark2 = np.mean(df_bi_101V["Current A"].iloc[200:330])
V101_light2 = np.mean(df_bi_101V["Current A"].iloc[365:496])
curr_bi_101.append(V101_light2-V101_dark2) 

V101_dark3 = np.mean(df_bi_101V["Current A"].iloc[520:710])
V101_light3 = np.mean(df_bi_101V["Current A"].iloc[740:850])
curr_bi_101.append(V101_light3-V101_dark3) 

V101_dark4 = np.mean(df_bi_101V["Current A"].iloc[900:1040])
V101_light4 = np.mean(df_bi_101V["Current A"].iloc[1100:1196])
curr_bi_101.append(V101_light4-V101_dark4)  

photo_current_bi = np.array([np.mean(curr_bi_21), np.mean(curr_bi_41), np.mean(curr_bi_81), np.mean(curr_bi_101)])
bias_bi_au = [21, 41, 81, 101]

plt.scatter(bias_bi_au, photo_current_bi*10**9)
plt.ylabel("Photocurrent [pA]")
plt.xlabel("Voltage [V]")

# fit 
popt_bi, pcov_bi = curve_fit(Hecht_eq_crystal, bias_bi_au, photo_current_bi, p0=[0.00000001, 122])
x_plot_bi = np.linspace(20,101, 1000)
plt.plot(x_plot_bi, Hecht_eq_crystal(x_plot_bi, *popt_bi)*10**9, color="coral")
plt.savefig("Hecht_fit_aubi.pdf")

print(f"mobility lifetime product is {popt_bi[0]}m^2/V with standard error {np.sqrt(pcov_bi[0][0])}m^2/V")

#molten film 

#df_mol_0V = pd.read_csv("Mol_filmi_0V_80uA.csv") # kinda bad
#df_mol_1V = pd.read_csv("Mol_filmi_1V_80uA.csv") # maybe lets not use this 
#df_ = pd.read_csv("Mol_filmi_2V_80uA.csv") # meh 
#df2 = pd.read_csv("Mol_filmi_3V_80uA.csv") # meh
 
df_mol_neg1 = pd.read_csv("Mol_filmi_neg1V_80uA.csv")
df_mol_neg3 = pd.read_csv("Mol_filmi_neg3V_80uA.csv")
df_mol_neg5 = pd.read_csv("Mol_filmi_neg5V_80uA.csv")
df_mol_neg7 = pd.read_csv("Mol_filmi_neg7V_80uA.csv")
df_mol_neg10 = pd.read_csv("Mol_filmi_neg10V_80uA.csv") 

df_mol_neg1["Current A"] = df_mol_neg1["Vsense 1 V (V)"]/R
df_mol_neg3["Current A"] = df_mol_neg3["Vsense 1 V (V)"]/R
df_mol_neg5["Current A"] = df_mol_neg5["Vsense 1 V (V)"]/R
df_mol_neg7["Current A"] = df_mol_neg7["Vsense 1 V (V)"]/R
df_mol_neg10["Current A"] = df_mol_neg10["Vsense 1 V (V)"]/R

plt.plot(df_mol_neg1["Current A"])
plt.title(r"Molten film -1 V 80 $\mu$A")
#plt.plot(df_mol_neg1["Current A"].iloc[660:735])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("molfilm_xrays_neg1V.pdf")
plt.savefig("molfilm_xrays_neg1V.eps")
plt.show() 

plt.plot(df_mol_neg3["Current A"])
plt.title(r"Molten film -3 V 80 $\mu$A")
#plt.plot(df_mol_neg3["Current A"].iloc[660:735])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("molfilm_xrays_neg3V.pdf")
plt.savefig("molfilm_xrays_neg3V.eps")
plt.show() 

plt.plot(df_mol_neg5["Current A"])
plt.title(r"Molten film -5 V 80 $\mu$A")
#plt.plot(df_mol_neg5["Current A"].iloc[660:735])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("molfilm_xrays_neg5V.pdf")
plt.savefig("molfilm_xrays_neg5V.eps")
plt.show() 

plt.plot(df_mol_neg7["Current A"])
plt.title(r"Molten film -7 V 80 $\mu$A")
#plt.plot(df_mol_neg7["Current A"].iloc[720:815])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("molfilm_xrays_neg7V.pdf")
plt.savefig("molfilm_xrays_neg7V.eps")
plt.show() 
 
plt.plot(df_mol_neg10["Current A"])
plt.title(r"Molten film -10 V 80 $\mu$A")
#plt.plot(df_mol_neg10["Current A"].iloc[630:730])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.savefig("molfilm_xrays_neg10V.pdf")
plt.savefig("molfilm_xrays_neg10V.eps")
plt.show() 

curr_mol_1 = []

neg1_dark1 = np.mean(df_mol_neg1["Current A"].iloc[0:82])
neg1_light1 = np.mean(df_mol_neg1["Current A"].iloc[90:145])
curr_mol_1.append(neg1_light1-neg1_dark1)

neg1_dark2 = np.mean(df_mol_neg1["Current A"].iloc[155:295])
neg1_light2 = np.mean(df_mol_neg1["Current A"].iloc[305:365])
curr_mol_1.append(neg1_light2-neg1_dark2)  

neg1_dark3 = np.mean(df_mol_neg1["Current A"].iloc[375:520])
neg1_light3 = np.mean(df_mol_neg1["Current A"].iloc[530:590])
curr_mol_1.append(neg1_light3-neg1_dark3)  

curr_mol_3 = []

neg3_dark1 = np.mean(df_mol_neg3["Current A"].iloc[0:45])
neg3_light1 = np.mean(df_mol_neg3["Current A"].iloc[70:150])
curr_mol_3.append(neg3_light1-neg3_dark1)

neg3_dark2 = np.mean(df_mol_neg3["Current A"].iloc[170:320])
neg3_light2 = np.mean(df_mol_neg3["Current A"].iloc[335:400])
curr_mol_3.append(neg3_light2-neg3_dark2)

neg3_dark3 = np.mean(df_mol_neg3["Current A"].iloc[420:560])
neg3_light3 = np.mean(df_mol_neg3["Current A"].iloc[570:690])
curr_mol_3.append(neg3_light3-neg3_dark3)

curr_mol_5 = []

neg5_dark1 = np.mean(df_mol_neg5["Current A"].iloc[0:45])
neg5_light1 = np.mean(df_mol_neg5["Current A"].iloc[60:179])
curr_mol_5.append(neg5_light1-neg5_dark1)

neg5_dark2 = np.mean(df_mol_neg5["Current A"].iloc[200:325])
neg5_light2 = np.mean(df_mol_neg5["Current A"].iloc[340:470])
curr_mol_5.append(neg5_light2-neg5_dark2)

neg5_dark3 = np.mean(df_mol_neg5["Current A"].iloc[500:630])
neg5_light3 = np.mean(df_mol_neg5["Current A"].iloc[660:735])
curr_mol_5.append(neg5_light3-neg5_dark3)

curr_mol_7 = []

neg7_dark1 = np.mean(df_mol_neg7["Current A"].iloc[0:50])
neg7_light1 = np.mean(df_mol_neg7["Current A"].iloc[80:170])
curr_mol_7.append(neg7_light1-neg7_dark1)

neg7_dark2 = np.mean(df_mol_neg7["Current A"].iloc[195:325])
neg7_light2 = np.mean(df_mol_neg7["Current A"].iloc[370:480])
curr_mol_7.append(neg7_light2-neg7_dark2)

neg7_dark3 = np.mean(df_mol_neg7["Current A"].iloc[520:665])
neg7_light3 = np.mean(df_mol_neg7["Current A"].iloc[720:815])
curr_mol_7.append(neg7_light3-neg7_dark3)

curr_mol_10 = []

neg10_dark1 = np.mean(df_mol_neg10["Current A"].iloc[0:45])
neg10_light1 = np.mean(df_mol_neg10["Current A"].iloc[75:160])
curr_mol_10.append(neg10_light1-neg10_dark1)

neg10_dark2 = np.mean(df_mol_neg10["Current A"].iloc[210:310])
neg10_light2 = np.mean(df_mol_neg10["Current A"].iloc[330:436])
curr_mol_10.append(neg10_light2-neg10_dark2)

neg10_dark3 = np.mean(df_mol_neg10["Current A"].iloc[470:600])
neg10_light3 = np.mean(df_mol_neg10["Current A"].iloc[630:730])
curr_mol_10.append(neg10_light3-neg10_dark3)

photo_current_mol = np.abs([np.mean(curr_mol_1), np.mean(curr_mol_3), np.mean(curr_mol_5), np.mean(curr_mol_7), np.mean(curr_mol_10)])
bias_mol = np.abs([-1, -3, -5, -7, -10])
plt.scatter(bias_mol, photo_current_mol*10**9)
plt.ylabel("Photocurrent [pA]")
plt.xlabel("Absolute Voltage [V]")

# fit 
popt_mol, pcov_mol = curve_fit(Hecht_eq_molfilm, bias_mol, photo_current_mol, p0=[0.0000002, 122])
x_plot_mol = np.linspace(0.1,10, 1000)
plt.plot(x_plot_mol, Hecht_eq_crystal(x_plot_mol, *popt_mol)*10**9, color="coral")
plt.savefig("Hecht_fit_molfilm.pdf")

print(f"mobility lifetime product is {popt_mol[0]}m^2/V with standard error {np.sqrt(pcov_mol[0][0])}m^2/V")

# Dose dependance for sensitivity measurement 

# area of contacts  
mol_film_area = 0.08**2 # pixel size cm 
half_contact_area = 0.35**2*np.pi/2 # size of Au Au / Au Bi contacts cm 

# at -3V
dose_au_20_uA = pd.read_csv("Au_Au_neg3V_20uA.csv")
dose_au_40_uA = pd.read_csv("Au_Au_neg3V_40uA.csv")
dose_au_60_uA = pd.read_csv("Au_Au_neg3V_60uA.csv")
dose_au_80_uA = pd.read_csv("Au_Au_neg3V_79uA_V2.csv")


dose_au_20_uA["Current A"] = dose_au_20_uA["Vsense 1 V (V)"]/R
dose_au_40_uA["Current A"] = dose_au_40_uA["Vsense 1 V (V)"]/R
dose_au_60_uA["Current A"] = dose_au_60_uA["Vsense 1 V (V)"]/R
dose_au_80_uA["Current A"] = dose_au_80_uA["Vsense 1 V (V)"]/R

plt.plot(dose_au_20_uA["Current A"])
plt.title(r"Dose gold gold -3 V 20 $\mu$A")
#plt.plot(dose_au_20_uA["Current A"].iloc[815:890])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

plt.plot(dose_au_40_uA["Current A"])
plt.title(r"Dose gold gold -3 V 40 $\mu$A")
#plt.plot(dose_au_40_uA["Current A"].iloc[880:977])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

plt.plot(dose_au_60_uA["Current A"])
plt.title(r"Dose gold gold -3 V 60 $\mu$A")
#plt.plot(dose_au_60_uA["Current A"].iloc[878:984])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

plt.plot(dose_au_80_uA["Current A"])
plt.title(r"Dose gold gold -3 V 80 $\mu$A")
#plt.plot(dose_au_80_uA["Current A"].iloc[75:168])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show() 

# lin reg function 

def lin_reg(x, k, b):
    return k*x+b

curr_dose_au_20 = []

dose_au_20_dark1 = np.mean(dose_au_20_uA["Current A"].iloc[0:60])
dose_au_20_light1 = np.mean(dose_au_20_uA["Current A"].iloc[70:160])
curr_dose_au_20.append(dose_au_20_light1-dose_au_20_dark1)

dose_au_20_dark2 = np.mean(dose_au_20_uA["Current A"].iloc[175:310])
dose_au_20_light2 = np.mean(dose_au_20_uA["Current A"].iloc[320:410])
curr_dose_au_20.append(dose_au_20_light2-dose_au_20_dark2)

dose_au_20_dark3 = np.mean(dose_au_20_uA["Current A"].iloc[420:580])
dose_au_20_light3 = np.mean(dose_au_20_uA["Current A"].iloc[594:690])
curr_dose_au_20.append(dose_au_20_light3-dose_au_20_dark3)

dose_au_20_dark4 = np.mean(dose_au_20_uA["Current A"].iloc[700:804])
dose_au_20_light4 = np.mean(dose_au_20_uA["Current A"].iloc[815:890])
curr_dose_au_20.append(dose_au_20_light4-dose_au_20_dark4)

#-----------------------------------------------------------------------

curr_dose_au_40 = []

dose_au_40_dark1 = np.mean(dose_au_40_uA["Current A"].iloc[0:77])
dose_au_40_light1 = np.mean(dose_au_40_uA["Current A"].iloc[94:200])
curr_dose_au_40.append(dose_au_40_light1-dose_au_40_dark1)

dose_au_40_dark2 = np.mean(dose_au_40_uA["Current A"].iloc[217:334])
dose_au_40_light2 = np.mean(dose_au_40_uA["Current A"].iloc[350:450])
curr_dose_au_40.append(dose_au_40_light2-dose_au_40_dark2)

dose_au_40_dark1 = np.mean(dose_au_40_uA["Current A"].iloc[470:606])
dose_au_40_light1 = np.mean(dose_au_40_uA["Current A"].iloc[627:726])
curr_dose_au_40.append(dose_au_40_light1-dose_au_40_dark1)

dose_au_40_dark1 = np.mean(dose_au_40_uA["Current A"].iloc[740:860])
dose_au_40_light1 = np.mean(dose_au_40_uA["Current A"].iloc[880:977])
curr_dose_au_40.append(dose_au_40_light1-dose_au_40_dark1)

#----------------------------------------------------------------------

curr_dose_au_60 = []

dose_au_60_dark1 = np.mean(dose_au_60_uA["Current A"].iloc[0:70])
dose_au_60_light1 = np.mean(dose_au_60_uA["Current A"].iloc[75:168])
curr_dose_au_60.append(dose_au_60_light1-dose_au_60_dark1)

dose_au_60_dark2 = np.mean(dose_au_60_uA["Current A"].iloc[198:320])
dose_au_60_light2 = np.mean(dose_au_60_uA["Current A"].iloc[340:420])
curr_dose_au_60.append(dose_au_60_light2-dose_au_60_dark2)

dose_au_60_dark3 = np.mean(dose_au_60_uA["Current A"].iloc[440:605])
dose_au_60_light3 = np.mean(dose_au_60_uA["Current A"].iloc[640:708])
curr_dose_au_60.append(dose_au_60_light3-dose_au_60_dark3)

dose_au_60_dark4 = np.mean(dose_au_60_uA["Current A"].iloc[725:859])
dose_au_60_light4 = np.mean(dose_au_60_uA["Current A"].iloc[878:984])
curr_dose_au_60.append(dose_au_60_light4-dose_au_60_dark4)

#---------------------------------------------------------------------

curr_dose_au_80 = []

dose_au_80_dark1 = np.mean(dose_au_80_uA["Current A"].iloc[0:60])
dose_au_80_light1 = np.mean(dose_au_80_uA["Current A"].iloc[80:168])
curr_dose_au_80.append(dose_au_80_light1-dose_au_60_dark1)

dose_au_80_dark2 = np.mean(dose_au_80_uA["Current A"].iloc[180:275])
dose_au_80_light2 = np.mean(dose_au_80_uA["Current A"].iloc[300:372])
curr_dose_au_80.append(dose_au_80_light2-dose_au_60_dark2)

dose_au_80_dark3 = np.mean(dose_au_80_uA["Current A"].iloc[390:500])
dose_au_80_light3 = np.mean(dose_au_80_uA["Current A"].iloc[520:610])
curr_dose_au_80.append(dose_au_80_light3-dose_au_60_dark3)

curr_dose_au_neg3V = np.array([np.mean(curr_dose_au_20), np.mean(curr_dose_au_40), np.mean(curr_dose_au_60), np.mean(curr_dose_au_80)])

curr_dens_au_neg3V = np.abs(curr_dose_au_neg3V/half_contact_area)

# linear regression to get the sensitivity  

dosage_au = np.array([12.48, 24.72, 36.84, 48.36])/60#Gy/s

plt.scatter(dosage_au, curr_dens_au_neg3V)
plt.xlabel(r"Dose rate [Gy s$^{-1}$]")
plt.ylabel(r"Current density [A cm$^{-2}$]")

#fit 
popt_dose_au, pcov_dose_au = curve_fit(lin_reg, dosage_au, curr_dens_au_neg3V) 
x_plot_dose_au = np.linspace(0, 1, 1000)
plt.plot(x_plot_dose_au, lin_reg(x_plot_dose_au, *popt_dose_au))
plt.savefig("dose_fit_auaa.pdf")

print(f"sensitivity is {popt_dose_au[0]} with stand error {np.sqrt(pcov_dose_au[0][0])} C/(Gy*cm^2)")

# gold bismuth

dose_bi_20_uA = pd.read_csv("Au_Bi_101V_20uA.csv")
dose_bi_40_uA = pd.read_csv("Au_Bi_101V_40uA.csv")
dose_bi_60_uA = pd.read_csv("Au_Bi_101V_60uA.csv")
dose_bi_80_uA = pd.read_csv("Au_Bi_101V_80uA.csv")

dose_bi_20_uA["Current A"] = dose_bi_20_uA["Vsense 1 V (V)"]/R
dose_bi_40_uA["Current A"] = dose_bi_40_uA["Vsense 1 V (V)"]/R
dose_bi_60_uA["Current A"] = dose_bi_60_uA["Vsense 1 V (V)"]/R
dose_bi_80_uA["Current A"] = dose_bi_80_uA["Vsense 1 V (V)"]/R

plt.plot(dose_bi_20_uA["Current A"])
plt.title(r"Dose gold bismuth 101 V 20 $\mu$A")
#plt.plot(dose_bi_20_uA["Current A"].iloc[470:558])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show() 

plt.plot(dose_bi_40_uA["Current A"])
plt.title(r"Dose gold bismuth 101 V 40 $\mu$A")
#plt.plot(dose_bi_40_uA["Current A"].iloc[535:625])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show() 

plt.plot(dose_bi_60_uA["Current A"])
plt.title(r"Dose gold bismuth 101 V 60 $\mu$A")
#plt.plot(dose_bi_60_uA["Current A"].iloc[601:740])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

plt.plot(dose_bi_80_uA["Current A"])
plt.title(r"Dose gold bismuth 101 V 80 $\mu$A")
#plt.plot(dose_bi_80_uA["Current A"].iloc[1100:1198])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

curr_dose_bi_20 = []

dose_bi_20_dark1 = np.mean(dose_bi_20_uA["Current A"].iloc[0:40])
dose_bi_20_light1 = np.mean(dose_bi_20_uA["Current A"].iloc[50:170])
curr_dose_bi_20.append(dose_bi_20_light1-dose_bi_20_dark1)

dose_bi_20_dark2 = np.mean(dose_bi_20_uA["Current A"].iloc[180:240])
dose_bi_20_light2 = np.mean(dose_bi_20_uA["Current A"].iloc[250:347])
curr_dose_bi_20.append(dose_bi_20_light2-dose_bi_20_dark2)

dose_bi_20_dark3 = np.mean(dose_bi_20_uA["Current A"].iloc[357:460])
dose_bi_20_light3 = np.mean(dose_bi_20_uA["Current A"].iloc[470:558])
curr_dose_bi_20.append(dose_bi_20_light3-dose_bi_20_dark3)

#----------------------------------------------------------------------

curr_dose_bi_40 = []

dose_bi_40_dark1 = np.mean(dose_bi_40_uA["Current A"].iloc[0:46])
dose_bi_40_light1 = np.mean(dose_bi_40_uA["Current A"].iloc[55:140])
curr_dose_bi_40.append(dose_bi_40_light1-dose_bi_40_dark1)

dose_bi_40_dark2 = np.mean(dose_bi_40_uA["Current A"].iloc[160:280])
dose_bi_40_light2 = np.mean(dose_bi_40_uA["Current A"].iloc[310:400])
curr_dose_bi_40.append(dose_bi_40_light2-dose_bi_40_dark2)

dose_bi_40_dark3 = np.mean(dose_bi_40_uA["Current A"].iloc[410:520])
dose_bi_40_light3 = np.mean(dose_bi_40_uA["Current A"].iloc[535:625])
curr_dose_bi_40.append(dose_bi_40_light3-dose_bi_40_dark3)

#-----------------------------------------------------------------------

curr_dose_bi_60 = []

dose_bi_60_dark1 = np.mean(dose_bi_60_uA["Current A"].iloc[0:60])
dose_bi_60_light1 = np.mean(dose_bi_60_uA["Current A"].iloc[70:189])
curr_dose_bi_60.append(dose_bi_60_light1-dose_bi_60_dark1)

dose_bi_60_dark2 = np.mean(dose_bi_60_uA["Current A"].iloc[200:310])
dose_bi_60_light2 = np.mean(dose_bi_60_uA["Current A"].iloc[330:439])
curr_dose_bi_60.append(dose_bi_60_light2-dose_bi_60_dark2)

dose_bi_60_dark3 = np.mean(dose_bi_60_uA["Current A"].iloc[460:590])
dose_bi_60_light3 = np.mean(dose_bi_60_uA["Current A"].iloc[601:740])
curr_dose_bi_60.append(dose_bi_60_light3-dose_bi_60_dark3)

#----------------------------------------------------------------------

curr_dose_bi_80 = []

dose_bi_80_dark1 = np.mean(dose_bi_80_uA["Current A"].iloc[150:330])
dose_bi_80_light1 = np.mean(dose_bi_80_uA["Current A"].iloc[370:490])
curr_dose_bi_80.append(dose_bi_80_light1-dose_bi_80_dark1)

dose_bi_80_dark2 = np.mean(dose_bi_80_uA["Current A"].iloc[528:700])
dose_bi_80_light2 = np.mean(dose_bi_80_uA["Current A"].iloc[745:840])
curr_dose_bi_80.append(dose_bi_80_light2-dose_bi_80_dark2)

dose_bi_80_dark3 = np.mean(dose_bi_80_uA["Current A"].iloc[900:1040])
dose_bi_80_light3 = np.mean(dose_bi_80_uA["Current A"].iloc[1100:1198])
curr_dose_bi_80.append(dose_bi_80_light3-dose_bi_80_dark3)

curr_dose_bi_101V = np.array([np.mean(curr_dose_bi_20), np.mean(curr_dose_bi_40), np.mean(curr_dose_bi_60), np.mean(curr_dose_bi_80)])

curr_dens_bi_101V = np.abs(curr_dose_bi_101V/half_contact_area)

# linear regression to get the sensitivity  

dosage_bi = np.array([12.48, 24.72, 36.84, 48.36])/60#Gy/min

plt.scatter(dosage_bi, curr_dens_bi_101V)
plt.xlabel(r"Dose rate [Gy s$^{-1}$]")
plt.ylabel(r"Current density [A cm$^{-2}$]")
#fit 
popt_dose_bi, pcov_dose_bi = curve_fit(lin_reg, dosage_bi, curr_dens_bi_101V) 
x_plot_dose_bi = np.linspace(0, 1, 1000)
plt.plot(x_plot_dose_bi, lin_reg(x_plot_dose_bi, *popt_dose_bi))

plt.savefig("dose_fit_aubi.pdf")
print(f"sensitivity is {popt_dose_bi[0]} with stand error {np.sqrt(pcov_dose_bi[0][0])}C/(Gy*cm^2)")

# molten film
dose_mol_20_uA = pd.read_csv("Mol_filmi_neg10V_20uA.csv")
dose_mol_40_uA = pd.read_csv("Mol_filmi_neg10V_40uA.csv")
dose_mol_60_uA = pd.read_csv("Mol_filmi_neg10V_60uA.csv")
dose_mol_80_uA = pd.read_csv("Mol_filmi_neg10V_80uA.csv")

dose_mol_20_uA["Current A"] = dose_mol_20_uA["Vsense 1 V (V)"]/R
dose_mol_40_uA["Current A"] = dose_mol_40_uA["Vsense 1 V (V)"]/R
dose_mol_60_uA["Current A"] = dose_mol_60_uA["Vsense 1 V (V)"]/R
dose_mol_80_uA["Current A"] = dose_mol_80_uA["Vsense 1 V (V)"]/R


plt.plot(dose_mol_20_uA["Current A"])
plt.title(r"Dose molten film -10 V 20 $\mu$A")
#plt.plot(dose_mol_20_uA["Current A"].iloc[1060:1190])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

plt.plot(dose_mol_40_uA["Current A"])
plt.title(r"Dose molten film -10 V 40 $\mu$A")
#plt.plot(dose_mol_40_uA["Current A"].iloc[690:810])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

plt.plot(dose_mol_60_uA["Current A"])
plt.title(r"Dose molten film -10 V 60 $\mu$A")
#plt.plot(dose_mol_60_uA["Current A"].iloc[680:798])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

plt.plot(dose_mol_80_uA["Current A"])
plt.title(r"Dose molten film -10 V 80 $\mu$A")
#plt.plot(dose_mol_80_uA["Current A"].iloc[615:730])
plt.xlabel("Number of measured points")
plt.ylabel("Current [A]")
plt.show()

curr_dose_mol_20 = []

dose_mol_20_dark1 = np.mean(dose_mol_20_uA["Current A"].iloc[0:56])
dose_mol_20_light1 = np.mean(dose_mol_20_uA["Current A"].iloc[70:255])
curr_dose_mol_20.append(dose_mol_20_light1-dose_mol_20_dark1)

dose_mol_20_dark2 = np.mean(dose_mol_20_uA["Current A"].iloc[280:430])
dose_mol_20_light2 = np.mean(dose_mol_20_uA["Current A"].iloc[450:534])
curr_dose_mol_20.append(dose_mol_20_light2-dose_mol_20_dark2)

dose_mol_20_dark3 = np.mean(dose_mol_20_uA["Current A"].iloc[560:670])
dose_mol_20_light3 = np.mean(dose_mol_20_uA["Current A"].iloc[700:830])
curr_dose_mol_20.append(dose_mol_20_light3-dose_mol_20_dark3)

dose_mol_20_dark4 = np.mean(dose_mol_20_uA["Current A"].iloc[850:1017])
dose_mol_20_light4 = np.mean(dose_mol_20_uA["Current A"].iloc[1060:1190])
curr_dose_mol_20.append(dose_mol_20_light4-dose_mol_20_dark4)

#-----------------------------------------------------------------------------

curr_dose_mol_40 = []


dose_mol_40_dark1 = np.mean(dose_mol_40_uA["Current A"].iloc[0:70])
dose_mol_40_light1 = np.mean(dose_mol_40_uA["Current A"].iloc[90:235])
curr_dose_mol_40.append(dose_mol_40_light1-dose_mol_40_dark1)

dose_mol_40_dark2 = np.mean(dose_mol_40_uA["Current A"].iloc[255:396])
dose_mol_40_light2 = np.mean(dose_mol_40_uA["Current A"].iloc[415:515])
curr_dose_mol_40.append(dose_mol_40_light2-dose_mol_40_dark2)

dose_mol_40_dark3 = np.mean(dose_mol_40_uA["Current A"].iloc[540:675])
dose_mol_40_light3 = np.mean(dose_mol_40_uA["Current A"].iloc[690:810])
curr_dose_mol_40.append(dose_mol_40_light3-dose_mol_40_dark3)

#-----------------------------------------------------------------------------

curr_dose_mol_60 = []


dose_mol_60_dark1 = np.mean(dose_mol_60_uA["Current A"].iloc[0:70])
dose_mol_60_light1 = np.mean(dose_mol_60_uA["Current A"].iloc[80:200])
curr_dose_mol_60.append(dose_mol_60_light1-dose_mol_60_dark1)

dose_mol_60_dark2 = np.mean(dose_mol_60_uA["Current A"].iloc[220:365])
dose_mol_60_light2 = np.mean(dose_mol_60_uA["Current A"].iloc[380:510])
curr_dose_mol_60.append(dose_mol_60_light2-dose_mol_60_dark2)

dose_mol_60_dark3 = np.mean(dose_mol_60_uA["Current A"].iloc[530:665])
dose_mol_60_light3 = np.mean(dose_mol_60_uA["Current A"].iloc[680:798])
curr_dose_mol_60.append(dose_mol_60_light3-dose_mol_60_dark3)

#------------------------------------------------------------------------------

curr_dose_mol_80 = []


dose_mol_80_dark1 = np.mean(dose_mol_80_uA["Current A"].iloc[0:45])
dose_mol_80_light1 = np.mean(dose_mol_80_uA["Current A"].iloc[60:160])
curr_dose_mol_80.append(dose_mol_80_light1-dose_mol_80_dark1)

dose_mol_80_dark2 = np.mean(dose_mol_80_uA["Current A"].iloc[180:310])
dose_mol_80_light2 = np.mean(dose_mol_80_uA["Current A"].iloc[335:440])
curr_dose_mol_80.append(dose_mol_80_light2-dose_mol_80_dark2)

dose_mol_80_dark3 = np.mean(dose_mol_80_uA["Current A"].iloc[460:600])
dose_mol_80_light3 = np.mean(dose_mol_80_uA["Current A"].iloc[615:730])
curr_dose_mol_80.append(dose_mol_80_light3-dose_mol_80_dark3)

curr_dose_mol_neg10V = np.array([np.mean(curr_dose_mol_20), np.mean(curr_dose_mol_40), np.mean(curr_dose_mol_60), np.mean(curr_dose_mol_80)])

curr_dens_mol_neg10V = np.abs(curr_dose_mol_neg10V/mol_film_area)

# linear regression to get the sensitivity  

dosage_mol = np.array([12.48, 24.72, 36.84, 48.36])/60# Gy/s

plt.scatter(dosage_mol, curr_dens_mol_neg10V)
plt.xlabel(r"Dose rate [Gy s$^{-1}$]")
plt.ylabel(r"Current density [A cm$^{-2}$]")
#fit 
popt_dose_mol, pcov_dose_mol = curve_fit(lin_reg, dosage_mol, curr_dens_mol_neg10V) 
x_plot_dose_mol = np.linspace(0, 1, 1000)
plt.plot(x_plot_dose_mol, lin_reg(x_plot_dose_mol, *popt_dose_mol))
plt.savefig("dose_fit_molfilm.pdf")

print(f"sensitivity is {popt_dose_mol[0]} with stand error {np.sqrt(pcov_dose_mol[0][0])}C/(Gy*cm^2)")