#program calling the class defined in the file stellar_core.py, StellarCore, and makes use of it in order to create and analyze data
#creates instances of the class modelling the density and temperature of the sun, as well as T=10^8 and checks accuracy 
#creates plots and prints data to the terminal analyzing energy output

from math import*
import numpy as np
import matplotlib.pyplot as plt
from stellar_core import StellarCore #importing the class StellarCore


T_8 = 10**8 ; T_1 = 1.57*10**7 #temperatures
rho = 1.62*10**5 #density constant of the sun
#definition of constants needed to perform calculations later:
k_B=1.3806*10**-23
h=6.6261*10**-34
e=1.6022*10**-19
eps=8.854*10**-12
u=1.660539*10**(-27)

#Sanity check for two temperatures, whereas T_1 simulates the temperature of the stellar core of the sun
t_c=[7.34*10**4, 1.10*10**0, 1.75*10**4, 1.23*10**-3, 4.35*10**-1, 1.27*10**5, 3.45*10**4] #set of expected values for the sanity check
err=0.01 #sanity check will pass if the relative error is smaller than this value
print("-------------------------------------------------------")
StellarCore(rho,T_8).sanity_check(t_c,err)
print("-------------------------------------------------------")
sun=StellarCore(rho,T_1)
t_c=[4.05*10**2,8.69*10**-9,4.87*10**-5,1.50*10**-6,5.30*10**-4,1.64*10**-6,9.18*10**-8] #set of expected values for the sanity check
sun.sanity_check(t_c,err)
print("-------------------------------------------------------")

sun=StellarCore(rho,T_1)
#energy lost from the different branches due to neutrinos is calculated and printed in the terminal
#calculated by Equations 6 and 1 from the report
total_energy=(4*sun.mH1-sun.mHe4) #total energy calculation
print("Branch:    Energy Lost to Neutrinos [%]:    Neutrino energy [MeV]:")
print(f"PP1:       {2*sun.E_PP0/(total_energy*sun.u_VJ)*100:.3f}                            {2*sun.E_PP0:.3f}")
print(f"PP2:       {(sun.E_PP0+sun.E_PP2)/(total_energy*sun.u_VJ)*100:.3f}                            {(sun.E_PP0+sun.E_PP2):.3f}")
print(f"PP3:       {(sun.E_PP0+sun.E_PP3)/(total_energy*sun.u_VJ)*100:.3f}                           {(sun.E_PP0+sun.E_PP3):.3f}")
print(f"CNO:       {(sun.E_CNO1+sun.E_CNO2)/(total_energy*sun.u_VJ)*100:.3f}                            {(sun.E_CNO1+sun.E_CNO2):.3f}")
print("-------------------------------------------------------")

#energy output per branch as a function of temperature, scaled and presented in a plot
colors = ["#88CEEB", "#CC6677", "#44AA99", "#DDCC77" , 'lightcoral', "#AA66BB", "gold", "forestgreen" , 'cadetblue', '#FF9966', "#661100", "gold","forestgreen"]
N=100001 #number of steps we wish to perform calculations for
#arrays to be containing energy output per branch of fusion-reactions
pp1=np.zeros(N)
pp2=np.zeros(N)
pp3=np.zeros(N)
cno=np.zeros(N)
T=np.linspace(10**4,10**9,N) #relevant intervall of temperatures
max_energy=[]
for i in range(len(pp1)):
    en=StellarCore(rho,T[i]) #an instance of the class for relevant density constant and temperature
    #calculations of energy outputs scaled by the total energy, saved in their respective arrays
    total_energy=en.Q_PP1+en.Q_PP2+en.Q_PP3+en.Q_CNO
    max_energy.append(total_energy)
    pp1[i]=en.Q_PP1/total_energy
    pp2[i]=en.Q_PP2/total_energy
    pp3[i]=en.Q_PP3/total_energy
    cno[i]=en.Q_CNO/total_energy 
#scaled energy plot
fig, axes = plt.subplots(1, 1, figsize=(10, 6))
axes.plot(T,pp1, color=colors[4], label="PP1", linewidth=1.6)
axes.plot(T,pp2, color=colors[3], label="PP3", linewidth=1.6)
axes.plot(T,pp3, color=colors[2], label="PP2", linewidth=1.6)
axes.plot(T,cno, color=colors[0], label="CNO", linewidth=1.6)
plt.xlabel(r"T [K]")
plt.ylabel(r"$\epsilon /\epsilon_{tot} \qquad [Jm^{-3}s^{-1}]$ ")
plt.xscale("log")
plt.legend()
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

#total energy plot
plt.plot(T,max_energy,color=colors[0],label="Total energy output")
plt.xlabel(r"T [K]")
plt.ylabel(r"$\epsilon_{tot} \qquad $ ")
plt.xscale("log")
plt.legend()
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

#Glamow peaks for the instance of a stellar core modelled with the temperature of the core of the sun
sun=StellarCore(rho,T_1)
N=5*10**4
E=np.linspace(10**-17,10**-13,N) #relevant energy interval
def peak_plot(E,T,Z1,Z2,m_1,m_2,lbl,clr): #function calculating the probability function, scaling by its max value, and plotting the energy against the function
    m=(m_1*m_2/(m_1+m_2))*u #\delta m 
    #peak defined and explained by the Equations 4 and 5 in the report
    peak=(np.exp(-E/(k_B*T)))*np.exp(-(np.sqrt(m/(2*E))*Z1*Z2*e**2*pi/(eps*h)))
    peak=peak/(np.sum(peak)*(abs((E[-1]-E[0])/len(E))))
    mx=np.max(peak)
    peak=peak/np.max(peak)
    def form(number):
        exp=int(np.floor(np.log10(abs(number))))
        coe=number/10**exp
        return fr"{coe:.3f}$\times 10^{exp}$"
    plt.plot(E, peak, label=f"{lbl}"+f",E={form(E[np.argmax(peak)])}J", linewidth=1.6,color=clr) #plots the scaled probability plot
#plots the scaled probability plot for each reaction of fusion occurring within the stellar core
peak_plot(E, T_1, 1, 1, sun.mH1, sun.mH1, "Eq 7", colors[0])
peak_plot(E, T_1, 1, 1, sun.mD2, sun.mH1, "Eq 8", colors[1])
peak_plot(E, T_1, 2, 2, sun.mHe3, sun.mHe3, "Eq 9", colors[2])
peak_plot(E, T_1, 2, 2, sun.mHe4, sun.mHe3, "Eq 10,13", colors[3])
peak_plot(E, T_1, 3, 1, sun.mLi7, sun.mH1, "Eq 12", colors[4])
peak_plot(E, T_1, 3, 1, sun.mBe7, sun.mH1, "Eq 14", colors[5])
peak_plot(E, T_1, 6, 1, sun.mC12, sun.mH1, "Eq 17", colors[6])
peak_plot(E, T_1, 6, 1, sun.mC13, sun.mH1, "Eq 19", colors[7])
peak_plot(E, T_1, 7, 1, sun.mN14, sun.mH1, "Eq 20", colors[8])
peak_plot(E, T_1, 7, 1, sun.mN15, sun.mH1, "Eq 22", colors[9])
plt.legend()
plt.xlabel(r"$E [J]$")
plt.ylabel(r"Proportionality ")
plt.xscale("log")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

