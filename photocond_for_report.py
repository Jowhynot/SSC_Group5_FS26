import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# variables needed for calculation 
area_sample = np.pi*0.35**2/2 # in cm^2 
area_ref = 0.785 # in cm^2 
C = 1.63*10**4 # in W^-1
R = 10 * 10**6 # 10 mega ohm 

def responsivity(sig_amp, sig_ph, sig_ref, area_sample=area_sample, C=C, R=R,area_ref=area_ref):
    return (C*sig_amp*area_ref*np.cos(sig_ph))/(R*sig_ref*area_sample)

def read_and_return_df(filename1, filename2):
    df_1 = pd.read_csv(filename1, sep="\t", names=["wavelength","idc" ,"phase", "reference", "idk2", "idk3", "idk4"])
    df_1.drop(columns=["idc", "idk2", "idk3", "idk4"], axis=1, inplace=True)
    df_1["phase"] = np.radians(df_1["phase"])
    df_2 = pd.read_csv(filename2, sep="\t", names=["signal", "mistake"], header=0).reset_index()
    df_1["signal"] = df_2["signal"]
    df = df_1.iloc[:-1]
    return df

#gold gold contact

df_au_au_0V_50nA = read_and_return_df("Au_Au_0V_50nA_ref_0dB_100mV", "Au_Au_0V_50nA_ref_0dB_100mV_1.txt")
df_au_au_0V_50nA["responsivity"] = np.abs(responsivity(sig_amp=df_au_au_0V_50nA.signal, sig_ph=df_au_au_0V_50nA.phase, sig_ref=df_au_au_0V_50nA.reference))
df_au_au_0V_50nA

df_au_au_neg2V_50nA = read_and_return_df("Au_Au_neg2V_50nA_ref_0dB_100mV", "Au_Au_neg2V_50nA_ref_0dB_100mV.txt")
df_au_au_neg2V_50nA["responsivity"] = np.abs(responsivity(sig_amp=df_au_au_neg2V_50nA.signal, sig_ph=df_au_au_neg2V_50nA.phase, sig_ref=df_au_au_neg2V_50nA.reference))
df_au_au_neg2V_50nA

df_au_au_neg5V_50nA = read_and_return_df("Au_Au_neg5V_50nA_ref_0dB_100mV", "Au_Au_neg5V_50nA_ref_0dB_100mV.txt")
df_au_au_neg5V_50nA["responsivity"] = np.abs(responsivity(sig_amp=df_au_au_neg5V_50nA.signal, sig_ph=df_au_au_neg5V_50nA.phase, sig_ref=df_au_au_neg5V_50nA.reference))
df_au_au_neg5V_50nA

df_au_au_neg8V_50nA = read_and_return_df("Au_Au_neg8V_100nA_ref_0dB_100mV", "Au_Au_neg8V_100nA_ref_0dB_100mV_1.txt")
df_au_au_neg8V_50nA["responsivity"] = np.abs(responsivity(sig_amp=df_au_au_neg8V_50nA.signal, sig_ph=df_au_au_neg8V_50nA.phase, sig_ref=df_au_au_neg8V_50nA.reference))
df_au_au_neg8V_50nA


plt.plot(df_au_au_0V_50nA.wavelength, df_au_au_0V_50nA.responsivity*10**9, label="0 V")
plt.plot(df_au_au_neg2V_50nA.wavelength, df_au_au_neg2V_50nA.responsivity*10**9, label="-2 V")
plt.plot(df_au_au_neg5V_50nA.wavelength, df_au_au_neg5V_50nA.responsivity*10**9, label="-5 V")
plt.plot(df_au_au_neg8V_50nA.wavelength, df_au_au_neg8V_50nA.responsivity*10**9, label="-8 V")
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"Responsivity/10$^9$ [A/W]")
plt.legend()
plt.show()

#bismuth gold 

df_au_bi_0V = read_and_return_df("Au_Bi_0V_50pA_ref_0dB_100mV", "Au_Bi_0V_50pA_ref_0dB_100mV_2.txt")
df_au_bi_0V["responsivity"] = np.abs(responsivity(sig_amp=df_au_bi_0V.signal, sig_ph=df_au_bi_0V.phase, sig_ref=df_au_bi_0V.reference))


df_au_bi_neg2V = read_and_return_df("Au_Bi_neg2V_100pA_ref_0dB_100mV", "Au_Bi_neg2V_100pA_ref_0dB_100mV_1.txt")
df_au_bi_neg2V["responsivity"] = np.abs(responsivity(sig_amp=df_au_bi_neg2V.signal, sig_ph=df_au_bi_neg2V.phase, sig_ref=df_au_bi_neg2V.reference))

df_au_bi_neg5V = read_and_return_df("Au_Bi_neg5V_500pA_ref_0dB_100mV", "Au_Bi_neg2V_100pA_ref_0dB_100mV_1.txt")
df_au_bi_neg5V["responsivity"] = np.abs(responsivity(sig_amp=df_au_bi_neg5V.signal, sig_ph=df_au_bi_neg5V.phase, sig_ref=df_au_bi_neg5V.reference))



df_au_bi_neg8V = read_and_return_df("Au_Bi_neg8V_500pA_ref_0dB_100mV","Au_Bi_neg8V_500pA_ref_0dB_100mV_1.txt")
df_au_bi_neg8V["responsivity"] = np.abs(responsivity(sig_amp=df_au_bi_neg8V.signal, sig_ph=df_au_bi_neg8V.phase, sig_ref=df_au_bi_neg8V.reference))


plt.plot(df_au_bi_0V.wavelength, df_au_bi_0V.responsivity*10**9, label="0 V")
plt.plot(df_au_bi_neg2V.wavelength, df_au_bi_neg2V.responsivity*10**9, label="-2 V")
plt.plot(df_au_bi_neg5V.wavelength, df_au_bi_neg5V.responsivity*10**9, label="-5 V")
plt.plot(df_au_bi_neg8V.wavelength, df_au_bi_neg8V.responsivity*10**9, label="-8 V")
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"Responsivity/10$^9$ [A/W]")
plt.legend()
plt.savefig("responsivity_photocond_au_bi.pdf")
plt.show()

#half moon

df_halfmoon_5V = read_and_return_df("hlafmoon_5V_50nA_ref_0dB_100mV","hlafmoon_5V_50nA_ref_0dB_100mV_1.txt")
df_halfmoon_5V["responsivity"] = np.abs(responsivity(sig_amp=df_halfmoon_5V.signal, sig_ph=df_halfmoon_5V.phase, sig_ref=df_halfmoon_5V.reference))

plt.plot(df_halfmoon_5V.wavelength, df_halfmoon_5V.responsivity*10**9, label="5 V")
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"Responsivity/10$^9$ [A/W]")
plt.legend()
plt.savefig("responsivity_photocond_halfmoon.pdf")
plt.show()

# molten film
df_molfilm_0V = read_and_return_df("Mol_film_0V_50pA_ref_0dB_100mV","Mol_film_0V_50pA_ref_0dB_100mV_1.txt")
df_molfilm_0V["responsivity"] = np.abs(responsivity(sig_amp=df_molfilm_0V.signal, sig_ph=df_molfilm_0V.phase, sig_ref=df_molfilm_0V.reference))

plt.plot(df_molfilm_0V.wavelength, df_molfilm_0V.responsivity*10**9, label="0 V")
plt.xlabel("Wavelength [nm]")
plt.ylabel(r"Responsivity/10$^9$ [A/W]")
plt.legend()
plt.savefig("responsivity_photocond_mol_film.pdf")
plt.show()
