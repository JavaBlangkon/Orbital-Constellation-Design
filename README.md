# Orbital-Constellation-Design
a thesis project to design low earth micro-satellites constellation for Air Traffic Monitoring system in Indonesia

SV2COE.py file consist code program to transform given state vector into classical orbital elements

COE2SV.py file consist code program to transform given classical orbital elements into state vector

twoBodyTool.py is a library tool consist the combination of function sv2coe, coe2sv, and derivative function F for the second ode integral function analysis

toolUtilizationExample.py is a code example of how to evaluate the integration of second ode problem as well as how to plot the relationships between orbital elements and the period

perturbationPlotTrial.py is the same method as toolUtilizationExample.py but now we are using the function of gravitational acceleration under perturbation due to J2

rateOfChangePerturbation.py is a method to evaluate the rate of change of each classical orbital elements using basic formula of mean orbital elements and compare it with the graph plotted of each classical orbital elements and analyze it using linear regression

coverageBelt.py is a method to calculate and estimating the coverage belt radius wideness of the satellite, the calculation depends strongly on orbital altitude and elevation angle

groundTrackPlot.py is the method to plot the ground track and orbital shape of the satellite, it includes the process to convert classical orbital elements into spherical coordinates, and defines the longitude and latitude of the satellite
