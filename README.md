# PMPS_PatankarSubroutine

[![DOI](https://zenodo.org/badge/259113602.svg)](https://zenodo.org/badge/latestdoi/259113602)

The subroutine used to solve the PMPRK scheme

In the directory of DataSamplingFullyDissipativeSystem, a fully dissipative/conservative system is solved.
In the directory of DataSamplingSemidissipativeSystem, all the equations except the last is dissipative/conservative.
The last equation representing the internal energy may not be dissipative, i.e., maybe increases due to the exothermal reactions.

Configurations at the top of "main.cpp" file:

int nSpecies = 100; #Number of species (plus an extra internale energy equation)

int nReactions = 100; #Number of reactions in the system

double deltatime = 1; #Default value, do not alter

double patankarCriticalValue = 1e-13; #The critical value used to stop the Newton iteration

int NTestedPerSample = 10; #How many random initial values generated for the Newton iteration used to test the uniqueness of                              #the solution

bool conservative = true;  #if true, the equations for the nSpecies+1 equations (DataSamplingFullyDissipativeSystem)
                           # or the nSpecies equations (DataSamplingSemiDissipativeSystem) is conservative.
                           
int maximumReactants = 10; #Maximum number of reactants considerred.

real AmplititudeOfInternalEnergyEquation = 100;  #The coefficients for the nSpecies density equations are drawn out from [0-1].
                                                 #The coefficients for the internale energy equation are drawn out from                                                          #[0,AmplititudeOfInternalEnergyEquation]
                                                 
int NTested = 10000;                             #Number of tested samples.
