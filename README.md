# Project 1 AST3310 : Model of a Stellar Engine 

This repository contains python code constructing a model for energy production within a stellar core, and this model is used to run simulations and exploring the model itself as well as the physics behind the energy production. The model for energy production is defined as a class StellarCore in the python file `stellar_core.py`, and the simulations are ran by calling and creating instances of this class in the python file `project1.py`.

-------------------------------------

### `stellar_core.py`: 
Contains a class `StellarCore` which takes a density constant $\rho$ and temperature T, and estimates the energy production at an arbitrary, fixed time. This energy production is a result of fusion of hydrogen into helium, through four different cycles of fusion reactions PP1, PP2, PP3 and CNO, derived and explained in the report. The class contains the following methods:

* `einsteins_principle(m_0,m_1):` which estimates energy from the einstein's mass-energy equivalency principle

* `number_density(mass_frac, Z):` estimates the number density of an arbitrary element Z

* `calculate_densities():` calls on the function `number_density` to calculate the number density of the elements relevant to the fusion reaction taking place within the core. The densities and isotopes are coupled and stored within a directory `n` so that the information can be extracted whenever

* `proportionality():` defines the relevant proportionality functions, for each of the relevant reactions of fusion. The proportionality is saved within a directory
  
* `reaction_rates(alter):` calculates the rates of reactions for of the relevant reactions of fusion from the calculated number densities and proportionality functions calculated and stored by `number_density()` and `proportionality()`. Takes 'alter' as an argument, which if alter=True is grounds for an alteration to the rate-calculation which takes into account realistic element production and consumption.

* `estimated_energy():` estimates the energy output of each reaction of fusion taking place within the stellar core from `einsteins_principle()`

* `cycle_relation():` calculates the total energy outputs from the reaction rates and energies calculated by `reaction_rates()` and `estimated_energy()`, and relates them to each of the cycles of fusion

* `sanity_check(t_c):` performs a test of the correctness of the model, for the temperatures $T=10^8$ or $T=1.57\times10^7$, and prints to the terminal if the sanity check is failed to be met 


### `project1.py`: 
Performs a number of simulations for different instances of the class StellarCore. First, a sanity check is performed for the density constant $\rho$ similar to that of the sun and the temperatures $T=10^8$ and $T=1.57\times10^7$. Furthermore, the energy lost from the different branches due to the neutrinos in the reactions is calculated in percent of the total energy output, and printed to the terminal. Lastly plots of enerrgy output as a function of temperatures, as well as plots of probability of certain energy plots - known as the Gamow plot - is created for each reaction of fusion taking place within the core is instituted, and saved in the seperate map 'figures' in form of a .pdf file. 
