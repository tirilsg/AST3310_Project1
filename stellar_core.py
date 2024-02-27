#program defining a class that models energy production within a stellar core, grounded in theory presented in the attatched report
#implements a number of relevant constants, as well as sizes like the mass of relevant isotopes created within the core
#defines a number of functions for number density of elements, reaction rates and energy output as a result of fusion
#a function testing the accuracy of the model is also implemented, sanity_check(), which can be called if needed

from math import*
import numpy as np
import matplotlib.pyplot as plt
import sys

class StellarCore:
    def __init__(self,rho,T):
        #definitions of relevant constants
        self.u_VJ=931.4943
        self.Na=6.022141*10**(23)
        self.g_si=1/(10**6*self.Na)
        self.MeV_J=10**6*1.602176*10**(-19)
        self.u=1.660539*10**(-27)
        self.rho=rho 
        self.T=T
        self.T9=T*10**(-9)

        #definition of relevant mass fractions
        self.X=0.7 ; self.Y3He=10**(-10) ; self.Y=0.29
        self.Z7Li=10**(-7) ; self.Z7Be=10**(-7) ; self.Z14N=10**(-11)

        #masses of isotopes relevant for the cycle reactions, units u
        self.mH1=1.007825 ; self.mD2=2.014102 ; self.mHe3=3.016029 ; self.mHe4=4.002603 
        self.mBe7=7.016929 ; self.mBe8=8.005305 ; self.mLi7=7.016004 ; self.mB8=8.024607 
        self.mC12=12.000000 ; self.mC13=13.003355 ; self.mN13=13.005739 ; self.mN14=14.003074
        self.mN15=15.000109 ; self.mO15=15.003066

        #energy lost to neutrinos per cycle, unit MeV
        self.E_PP0=0.265 ; self.E_PP2=0.815 ; self.E_PP3=6.711
        self.E_CNO1=0.707 ; self.E_CNO2=0.997

        #calculations of the energy output for each reaction cycle
        self.calculate_densities()
        self.proportionality()
        self.reaction_rates(True)
        self.estimated_energy()
        self.cycle_relation()

    def einsteins_principle(self,m_0,m_1): #Einstein's mass-energy equivalency principle, Equation 1 in the report
        return (m_0-m_1)*self.u_VJ

    def number_density(self, mass_frac, Z): #number density of an arbitrary element Z, Equation 2 in the report
        return mass_frac*self.rho/(Z*self.u)
    
    def calculate_densities(self): #calculation of the number density for all relevant isotopes
        self.n={} #definition of a dictionary that will be used to couple and store densities and isotopes
        self.n["1H"]=self.number_density(self.X,1)
        self.n["3He"]=self.number_density(self.Y3He,3)
        self.n["4He"]=self.number_density(self.Y,4)
        self.n["7Li"]=self.number_density(self.Z7Li,7)
        self.n["7Be"]=self.number_density(self.Z7Be,7)
        self.n["14N"]=self.number_density(self.Z14N,14)
        self.n["e"]=(self.n["1H"]+2*self.n["3He"]+2*self.n["4He"])+(3*self.n["7Li"]+3*self.n["7Be"]+7*self.n["14N"])
        return self.n
    
    def proportionality(self): #definition of all relevant proportionality functions, presented in the Table III in the report
        self.NA_lambda={} #dictionary that will contain data for proportionality coupled with the relevant fusion-reaction "name"
        #definition of constants needed in the estimation of proportionality functions
        T9=self.T9 
        T9S=T9/(1+4.95*10**(-2)*T9)
        T9SS=T9/(1+0.759*T9)
        self.NA_lambda["lpp"]=self.g_si*4.01*10**(-15)*T9**(-2/3)*np.exp(-3.380*T9**(-1/3))*(1+0.123*T9**(1/3)+1.09*T9**(2/3)+0.938*T9)
        self.NA_lambda["l33"]=self.g_si*6.04*10**10*T9**(-2/3)*np.exp(-12.276*T9**(-1/3))*(1+0.034*T9**(1/3)-0.522*T9**(2/3)-0.124*T9+0.353*T9**(4/3)+0.213*T9**(5/3))
        self.NA_lambda["l34"]=self.g_si*5.61*10**6*T9S**(5/6)*T9**(-3/2)*np.exp(-12.826*T9S**(-1/3))
        lambda_e7=(1.34*10**(-10)*T9**(-1/2)*(1-0.537*T9**(1/3)+3.86*T9**(2/3)+0.0027*T9**(-1)*np.exp(2.515*10**(-3)*T9**(-1))))*self.g_si
        if self.T<10**6: #test that insures the existence of the upper limit for Li_3^7 is considered
            if lambda_e7 > (1.57*10**(-7))*10**-6/(self.n["e"]*self.Na):
                lambda_e7 = (1.57*10**(-7))*10**-6/(self.n["e"]*self.Na)
        self.NA_lambda["le7"]=lambda_e7
        self.NA_lambda["l17l"]=self.g_si*(1.096*10**9*T9**(-2/3)*np.exp(-8.472*T9**(-1/3))-4.830*10**8*T9SS**(5/6)*T9**(-3/2)*np.exp(-8.472*T9SS**(-1/3))+1.06*10**10*T9**(-2/3)*np.exp(-30.442*T9**(-1)))
        self.NA_lambda["l17"]=self.g_si*(3.11*10**5*T9**(-2/3)*np.exp(-10.262*T9**(-1/3))+2.53*10**3*T9**(-3/2)*np.exp(-7.306*T9**(-1)))
        self.NA_lambda["lp14"]=self.g_si*((4.90*10**7*T9**(-2/3)*np.exp(-15.228*T9**(-1/3)-0.092*T9**2)*(1+0.027*T9**(1/3)-0.778*T9**(2/3)-0.149*T9+0.261*T9**(4/3)+0.127*T9**(5/3)))+2.37*10**3*T9**(-3/2)*np.exp(-3.011*T9**(-1))+2.19*10**4*np.exp(-12.53*T9**(-1)))
        return self.NA_lambda

    def reaction_rates(self, alter): #definition of rate of reactions, presented and explained as Equation 3 in the report
        self.rates={} #dictionary to contain reaction rates of each relevant fusion-reaction
        self.rates["lpp"] = self.n["1H"]*self.n["1H"]/(self.rho*(1+1))*self.NA_lambda["lpp"]
        self.rates["l33"] = self.n["3He"]*self.n["3He"]/(self.rho*(1+1))*self.NA_lambda["l33"]
        self.rates["l34"] = self.n["3He"]*self.n["4He"]/(self.rho)*self.NA_lambda["l34"]
        self.rates["le7"] = self.n["e"]*self.n["7Be"]/(self.rho)*self.NA_lambda["le7"]
        self.rates["l17l"] = self.n["7Li"]*self.n["1H"]/(self.rho)*self.NA_lambda["l17l"]
        self.rates["l17"] = self.n["7Be"]*self.n["1H"]/(self.rho)*self.NA_lambda["l17"]
        self.rates["lp14"] = self.n["14N"]*self.n["1H"]/(self.rho)*self.NA_lambda["lp14"]
        #if the argument alter is True, reaction rates are altered to make sure energy consumption in the stellar core
        #does not exceed the rate of production of elements, in accordance with expressions presented in Table IV in the report
        if alter==True: 
            const=self.rates["lpp"]/(2*self.rates["l33"]+self.rates["l34"])
            #if the rate of production of 2*l33 (lpp needs to happen 2 times in order for l33 to happen) and l34 is higher than lpp, rates needs to be adjusted
            if (2*self.rates["l33"]+self.rates["l34"]) > self.rates["lpp"]:
                self.rates["l33"] = const*self.rates["l33"]
                self.rates["l34"] = const*self.rates["l34"]
            const=self.rates["l34"]/(self.rates["le7"]+self.rates["l17"])
            #the fallout of l34 can be either reaction l17 and le7, and if the sum of these consumption rates exceeds the rate of production l34;
            if self.rates["le7"]+self.rates["l17"] > self.rates["l34"]:
                self.rates["le7"] = const*self.rates["le7"]
                self.rates["l17"] = const*self.rates["l17"]
            #if more of l17l happens more frequently than le7, the rates need to be adjusted to be realistic
            if self.rates["l17l"]>self.rates["le7"] :
                self.rates["l17l"] = self.rates["le7"]
        return self.rates

    def estimated_energy(self): #definition of total energy output of each reaction of fusion
        self.E={}
        self.E["l1"] = self.MeV_J*(self.einsteins_principle(2*self.mH1 , self.mD2) - self.E_PP0)
        self.E["l1d"] = self.MeV_J*self.einsteins_principle(self.mD2 + self.mH1 , self.mHe3)
        self.E["l33"] = self.MeV_J*self.einsteins_principle(2*self.mHe3 , self.mHe4 + 2*self.mH1)
        self.E["l34"] = self.MeV_J*self.einsteins_principle(self.mHe3 + self.mHe4 , self.mBe7)
        self.E["le7"] = self.MeV_J*(self.einsteins_principle(self.mBe7 , self.mLi7) - self.E_PP2)
        self.E["l17l"] = self.MeV_J*self.einsteins_principle(self.mLi7 + self.mH1 , 2*self.mHe4)
        self.E["l17"] = self.MeV_J*self.einsteins_principle(self.mBe7 + self.mH1 , self.mB8)
        self.E["dc"] = self.MeV_J*(self.einsteins_principle(self.mB8 , 2*self.mHe4) - self.E_PP3)
        self.E["lp12"] = self.MeV_J*self.einsteins_principle(self.mC12 + self.mH1 , self.mN13)
        self.E["l13"] = self.MeV_J*(self.einsteins_principle(self.mN13 , self.mC13) - self.E_CNO1)
        self.E["lp13"] = self.MeV_J*self.einsteins_principle(self.mC13 + self.mH1 , self.mN14)
        self.E["lp14"] = self.MeV_J*self.einsteins_principle(self.mN14 + self.mH1 , self.mO15)
        self.E["l15"] = self.MeV_J*(self.einsteins_principle(self.mO15 , self.mN15) - self.E_CNO2)
        self.E["lp15"] = self.MeV_J*self.einsteins_principle(self.mN15 + self.mH1 , self.mC12 + self.mHe4)
        return self.E
    
    def cycle_relation(self): #relation between each fusion-cycle and its energy output calculated by Equation 6 in the report
        self.Q_PP1 = self.rates["l33"] * (self.E["l33"] + 2*(self.E["l1"] + self.E["l1d"]))
        self.Q_PP2 = self.rates["l34"] * (self.E["l34"] + (self.E["l1"] +self.E["l1d"])) + self.rates["le7"]*(self.E["le7"] + self.E["l17l"])
        self.Q_PP3 = self.rates["l17"] * (self.E["l17"]+self.E["dc"])+self.rates["l34"]*(self.E["l34"]+(self.E["l1"] + self.E["l1d"]))
        self.Q_CNO = self.rates["lp14"] * (self.E["lp12"]+self.E["l13"]+self.E["lp13"]+self.E["lp14"]+self.E["lp15"]+self.E["l15"])
        return self.Q_PP1, self.Q_PP2, self.Q_PP3, self.Q_CNO
    
    def sanity_check(self,t_v,err): #a function that will check the accuracy of the developed model
        print(f"Sanity check for temperature T={self.T} and density rho={self.rho}")
        print("Expected value:   Estimated value:   Relative Error:")
        e_v=np.array([
            self.rates["lpp"]*(self.E["l1"] + self.E["l1d"])*self.rho,
            self.rates["l33"]*self.E["l33"]*self.rho,
            self.rates["l34"]*self.E["l34"]*self.rho,
            self.rates["le7"]*self.E["le7"]*self.rho,
            self.rates["l17l"]*self.E["l17l"]*self.rho,
            self.rates["l17"]*(self.E["l17"]+self.E["dc"])*self.rho,
            self.rates["lp14"]*(self.E["lp12"]+self.E["l13"]+self.E["lp13"]+self.E["lp14"]+self.E["lp15"]+self.E["l15"])*self.rho])
        #prints the expected and estimated value for the energy production per unit mass, as well as the relative error is printed to the terminal
        list=[]
        for i in range(len(t_v)):
            rel=abs((e_v[i]-t_v[i])/t_v[i])
            if rel>err:
                list.append(i)
            print(f"{t_v[i]:.3}          {e_v[i]:.3}              {rel:.3}")
        if len(list)!=0:
            for i in range(len(list)):
                print(f"The sanity check is failed for Sanity Check Equation {list[i]+1}")
            #sys.exit() #exits entire program, meaning plots and other code trying to be ran after calling a sanity check that fails will not be able to run


