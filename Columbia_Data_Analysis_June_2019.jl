#Columbia_Data_Analysis_June_2019
#Thomas Moore, June 2019
#Run with Julia v1.0.2

#---------------------------------------------------------#
#Introduction to File
#---------------------------------------------------------#

#In this file, we will analyse some of the data sent by Guanhe
#on the absorption of CO2 at 1 bar into SIPs containing various
#weight fractions of NOHM-PEI (NPEI).

#We have absorption data for 5 different SIPs, containing 9.8wt%,
#29.4wt%, 33.3wt%, 40wt% and 58.8wt% NPEI inside Tegorad 2650 PDMS.

#The goal of this document is to analyse the data, in particular
#to determine the consistency of the *equilibrium* absorption
#capacities for the different materials, and to explore various
#approaches to modelling the *kinetics* of the CO2 uptake. The
#kinetics of the uptake are complicated by the presence of an
#autocatalytic-type reaction, which causes the system to be strongly
#*reaction-controlled* at the *start* of the absorption process. We
#will consider various approaches to modelling this behaviour.

#We will begin by calling relevant packages, and we will then make
#plots of the raw data from Guanhe. We will then refine this data so
#it only contains the absorption curve we are interested in, and we
#will make plots of these curves. After this, we will conduct some
#calculations on the equilibrium capacity and the kinetic uptake
#of CO2 into the material.

#---------------------------------------------------------#
#Call Relevant Packages
#---------------------------------------------------------#

using Plots
using Plots.PlotMeasures
using DelimitedFiles
using LaTeXStrings
using Statistics
using LsqFit

#---------------------------------------------------------#
#Load Absorption Data
#---------------------------------------------------------#

#The absorption data is contained in the CO2_Uptake_Data.txt file,
#which was previously created from the excel data sent by Guanhe.
#In this file, the first and second columns are time data and
#SIP mass data for the 9.8wt% case, the third and fourth colums are
#time and SIP mass data for the third and fourth columns, etc.

#The units are:
# - Time in Hours
# - Mass in Grams

#It may be necessary to run something like:
#cd("/Users/thomasmoore/Dropbox/Papers/NOHMandSIPs/NOHMsgit/")
RawUptakeData = readdlm("CO2_Uptake_Data.txt")

#Create Plot of Raw Absorption Data Curves
Plots.scalefontsizes(); Plots.scalefontsizes(1.4)
plot(RawUptakeData[:,1],RawUptakeData[:,2], color="black", linestyle = :solid, xlabel = L"\mathrm{Time\ (hr)}", ylabel = L"\mathrm{Mass\ (g)}", title = L"\mathrm{Raw\ CO_2\ Uptake\ Data}", label = "10 wt%")
plot!(RawUptakeData[:,3],RawUptakeData[:,4], color="blue", linestyle = :solid, label = "29 wt%")
plot!(RawUptakeData[:,5],RawUptakeData[:,6], color="red", linestyle = :solid, label = "33 wt%")
plot!(RawUptakeData[:,7],RawUptakeData[:,8], color="green", linestyle = :solid, label = "49 wt%")
plot!(RawUptakeData[:,9],RawUptakeData[:,10], color="grey", linestyle = :solid, label = "59 wt%")

#Create array containing just the absorption isotherms. This array
#will have the same structure as the RawUptakeData array: the first
#and second columns will contain time and mass data for the 10wt%
#curves, the third and fourth columns will contain time and mass
#data for the 30wt% curve, etc.

#In each case, CO2 uptake begins at a time step just before 3 hrs,
#at which point the temperature, which has been cooling since the
#2 hr mark, has just begun to stabilise. The exact step at which
#CO2 absorption begins is row 3349 in RawUptakeData, and continues
#for exactly 5 hours to row 8974. All timesteps in RawUptakeData
#are the same (3.2 seconds, or 1125 samples per hour), and all
#curves will be shifted to start at t = 0 at the begining of absorption.
#Likewise, all absorption isotherms will be shifted so that they
#begin with initial mass of 0.

#We will also convert the units to:

# - Time in seconds
# - Mass in kilograms.

SorptionIsothermData = RawUptakeData[3349:8974,:]
for i in 1:2:9 SorptionIsothermData[:,i] = (SorptionIsothermData[:,i] .- SorptionIsothermData[1,i]) .* 3600 end
for i in 2:2:10 SorptionIsothermData[:,i] = (SorptionIsothermData[:,i] .- SorptionIsothermData[1,i]) ./ 1000  end

#Create Plot of Absorption Isotherms
Plots.scalefontsizes(); Plots.scalefontsizes(1.4)
plot(SorptionIsothermData[:,1],SorptionIsothermData[:,2], color="black", linestyle = :solid, xlabel = L"\mathrm{Time\ (s)}", ylabel = L"\mathrm{Mass\ Change\ (kg)}", title = L"\mathrm{CO_2\ Absorption\ Isotherms}", label = "10 wt%", legend = :topleft, xlims = (0,18000), ylims = (0.0,1.0e-3))
plot!(SorptionIsothermData[:,3],SorptionIsothermData[:,4], color="blue", linestyle = :solid, label = "29 wt%")
plot!(SorptionIsothermData[:,5],SorptionIsothermData[:,6], color="red", linestyle = :solid, label = "33 wt%")
plot!(SorptionIsothermData[:,7],SorptionIsothermData[:,8], color="green", linestyle = :solid, label = "49 wt%")
plot!(SorptionIsothermData[:,9],SorptionIsothermData[:,10], color="grey", linestyle = :solid, label = "59 wt%")

#Plots.scalefontsizes(); Plots.scalefontsizes(1.4)
#plot(sqrt.(SorptionIsothermData[:,1]),SorptionIsothermData[:,2], color="black", linestyle = :solid, xlabel = L"\mathrm{Time\ (s)}", ylabel = L"\mathrm{Mass\ Change\ (kg)}", title = L"\mathrm{CO_2\ Absorption\ Isotherms}", label = "10 wt%", legend = :topleft, xlims = (0,sqrt(18000)), ylims = (0.0,1.0e-3))
#plot!(sqrt.(SorptionIsothermData[:,3]),SorptionIsothermData[:,4], color="blue", linestyle = :solid, label = "29 wt%")
#plot!(sqrt.(SorptionIsothermData[:,5]),SorptionIsothermData[:,6], color="red", linestyle = :solid, label = "33 wt%")
#plot!(sqrt.(SorptionIsothermData[:,7]),SorptionIsothermData[:,8], color="green", linestyle = :solid, label = "49 wt%")
#plot!(sqrt.(SorptionIsothermData[:,9]),SorptionIsothermData[:,10], color="grey", linestyle = :solid, label = "59 wt%")

#---------------------------------------------------------#
#Calculate and Compare Theoretical Capacities
#---------------------------------------------------------#

#We will calculate the amount of CO2 stored in the NPEI phase of
#each of these materials. The mass of each film is given by the mass
#measured in the TGA at the very beginning of absorption, after all
#water has evaporated off; this will be the mass in row 3349 in
#RawUptakeData. We will define the NPEI mass fractions in MassFractions,
#and we will define the film masses in kg in SIPMasses. We will then calculate
#the moles of CO2 absorbed per gram of NPEI from the beginning to
#the end of the 5 hour absorption run. Visually, it looks like this will
#be the saturation concentration within the SIP containing 10wt% NPEI,
#but the other SIPs don't quite looke like they've reached saturation.

#In order to account for the CO2 solubtility in the PDMS, we will use
#a solubility for pure CO2 in PDMS at 25C, interpolated from Scholes et al.
#(2010) from the Journal of Membrane Science, Table 3. This gives
#S = 5.94e-4 mol/m3.Pa. We will assume a density of PDMS of 965 kg/m3.

#We store the NPEI mass fractions in MassFractions
MassFractions = [0.098, 0.294, 0.333, 0.49, 0.588]

#We store the SIP film masses in SIPMasses
SIPMasses = RawUptakeData[3349,2:2:10]*1e-3                  #kg

#Define density and CO2 solubtility within PDMS.
SPDMS = 5.94e-4                                             #mol/m3.Pa
ρPDMS = 965                                                 #kg/m3

#Define CO2 Pressure
PCO2 = 1e5                                                  #Pa

#Calculate moles of CO2 stored in 1kg of NPEI after absorption for each of
#the 5 absorption runs.
MaxMolesCO2PerKgNPEI = [((SorptionIsothermData[end,2*i] - SorptionIsothermData[1,2*i])/0.044 - SIPMasses[i]*(1-MassFractions[i])/ρPDMS*SPDMS*PCO2)/(SIPMasses[i]*MassFractions[i]) for i = 1:5]

#Each of the various absorption isotherms predicts an `equilibrium` (after 5
#hours) CO2 capacity of between ~5-6. The final, 58.8wt% sample gives a
#slightly lower capacity, however this is unsurprising given that the isotherm
#looks as thought it has not reached complete equilibrium. The same could be
#said for the other curves to a lesser extent.

#We see the relatively consistent results from the following bar chart.
bar(["10","29","33","49","59"],MaxMolesCO2PerKgNPEI,legend=false, xlabel="NOHM wt%",ylabel="mol CO2/kg NOHM")

#Ignoring the 59wt% run, the mean capacity appears to be around 5.5mol/kg.
mean(MaxMolesCO2PerKgNPEI[1:end-1])

#---------------------------------------------------------#
#Analysis of Kinetics; Naive Approach.
#---------------------------------------------------------#

#To begin, we will consider a naive approach, in which we don't account for
#autocatalytic effects. We will use the model developed in Section 4.8 of
#my thesis. This will require:
# - CO2 Partial Pressure, pCO2 = 1 bar = 100,000 Pa
# - Permeability of CO2 = 3250 barrer = 1.1e-12 mol/m.Pa.s
# - Thickness of SIP, L = 1mm = 1e-3m
# - Solubility within PDMS = 5.94e-4 mol/Pa.m3
# - Capacity within NPEI = 5.5 mol/kg
# - We set the density of NPEI at 1250 kg/m3, however this cancels during the calculation, so the exact value is irrelevant.

pCO2 = 1e5                                                      #Pressure of CO2 Pa
PPDMS = 1.1e-12                                                 #CO2 Permeability of PDMS mol/m.Pa.s
PNPEI = 0                                                       #CO2 Permeability of NPEI mol/m.Pa.s
L = 1.9e-3                                                        #Sheet Thickness m
CapNPEI = 5.5                                                   #CO2 Capacity of NPEI mol/kg
ρNPEI = 1250                                                    #Density of NPEI kg/m3
ρPDMS = 965                                                     #Density of PDMS kg/m3

#We *could* calculate C1 using the equation from Section 4.7 of the thesis:
εMass =  1-MassFractions[1]                                     #Mass Fraction of PDMS
ε = εMass/ρPDMS / (εMass/ρPDMS + (1-εMass)/ρNPEI)
N = ε*SPDMS*pCO2 + (1-ε)*ρNPEI*CapNPEI
C1 = 2*pCO2/(N*L^2)*((1-ε)*PNPEI + ε*PPDMS)

#However, we assume PNPEI = 0 and combine these into one equation, which we
#calculate for each material:
C1_vals = [2*pCO2/L^2*PPDMS * (SPDMS*pCO2 + MassFractions[i]/(1-MassFractions[i])*ρPDMS*CapNPEI)^-1 for i = 1:5]

#In order to compare the model with the experimental CO2 absorption isotherm data,
#we need a theoretical estimate for the material capacity. We calculate
#the capacity in moles per kg.
SIPCapacity = [MassFractions[i]*CapNPEI + (1-MassFractions[i])/ρPDMS*SPDMS*pCO2 for i in 1:5]        #mol/kg

#Creaction functions to calculate predicted uptake curves
nCO2_10(t) = (t < 1/C1_vals[1] ? sqrt(C1_vals[1]*t) : 1)*SIPMasses[1]*SIPCapacity[1]
nCO2_29(t) = (t < 1/C1_vals[2] ? sqrt(C1_vals[2]*t) : 1)*SIPMasses[2]*SIPCapacity[2]
nCO2_33(t) = (t < 1/C1_vals[3] ? sqrt(C1_vals[3]*t) : 1)*SIPMasses[3]*SIPCapacity[3]
nCO2_49(t) = (t < 1/C1_vals[4] ? sqrt(C1_vals[4]*t) : 1)*SIPMasses[4]*SIPCapacity[4]
nCO2_58(t) = (t < 1/C1_vals[5] ? sqrt(C1_vals[5]*t) : 1)*SIPMasses[5]*SIPCapacity[5]

#Create plots comparing model with CO2 Absorption Isotherms
Plots.scalefontsizes(); Plots.scalefontsizes(1)
plot(SorptionIsothermData[:,1],SorptionIsothermData[:,2]/.044, color="black", linestyle = :solid, xlabel = L"\mathrm{Time\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ Isotherms}", ylims = (0,0.025), xlims = (0,18000),label = "10 wt%", legend = :topleft)
plot!(SorptionIsothermData[:,1],nCO2_10.(SorptionIsothermData[:,1]), color="black", linestyle = :dash, label = "10 wt% Model")
plot!(SorptionIsothermData[:,3],SorptionIsothermData[:,4]/0.044, color="blue", linestyle = :solid, label = "29 wt%")
plot!(SorptionIsothermData[:,3],nCO2_29.(SorptionIsothermData[:,3]), color="blue", linestyle = :dash, label = "29 wt% Model")
plot!(SorptionIsothermData[:,5],SorptionIsothermData[:,6]/0.044, color="red", linestyle = :solid, label = "33 wt%")
plot!(SorptionIsothermData[:,5],nCO2_33.(SorptionIsothermData[:,5]), color="red", linestyle = :dash, label = "33 wt% Model")
plot!(SorptionIsothermData[:,7],SorptionIsothermData[:,8]/0.044, color="green", linestyle = :solid, label = "49 wt%")
plot!(SorptionIsothermData[:,7],nCO2_49.(SorptionIsothermData[:,7]), color="green", linestyle = :dash, label = "49 wt% Model")
plot!(SorptionIsothermData[:,9],SorptionIsothermData[:,10]/0.044, color="grey", linestyle = :solid, label = "59 wt%")
plot!(SorptionIsothermData[:,9],nCO2_58.(SorptionIsothermData[:,9]), color="grey", linestyle = :dash, label = "59 wt% Model")

#---------------------------------------------------------#
#Plot of c vs sqrt(t)
#---------------------------------------------------------#

#We have been assuming to this point that after the small period of time in which
#we have autocatalytic behaviour, the absorption is diffuction controlled.
#If this is indeed the case, then we would expect a c vs sqrt(t) or c vs sqrt(t-t0)
#plot to be linear for much of the absorption regime. We will now create
#these plots to determine whether absorption is indeed ever diffusion controlled.

#Create Plots of c vs sqrt(t)
Plots.scalefontsizes(); Plots.scalefontsizes(1)
plot(sqrt.(SorptionIsothermData[:,1]),SorptionIsothermData[:,2]/.044, color="black", linestyle = :solid, xlabel = L"\sqrt{t}\mathrm{\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t}", ylims = (0,0.025), xlims = (0,140),label = "10 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[:,3]),SorptionIsothermData[:,4]/.044, color="blue", linestyle = :solid, label = "29 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[:,5]),SorptionIsothermData[:,6]/.044, color="red", linestyle = :solid, label = "33 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[:,7]),SorptionIsothermData[:,8]/.044, color="green", linestyle = :solid, label = "49 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[:,9]),SorptionIsothermData[:,10]/.044, color="grey", linestyle = :solid, label = "58 wt%", legend = :topleft)

#We now shift each of these in time, so that the straight line passes through the origin.
Plots.scalefontsizes(); Plots.scalefontsizes(1)
plot(sqrt.(SorptionIsothermData[1:end-71,1]),SorptionIsothermData[72:end,2]/.044, color="black", linestyle = :solid, xlabel = L"\sqrt{t}\mathrm{\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t}", ylims = (0,0.025), xlims = (0,140),label = "10 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[1:end-319,3]),SorptionIsothermData[320:end,4]/.044, color="blue", linestyle = :solid, label = "29 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[1:end-319,5]),SorptionIsothermData[320:end,6]/.044, color="red", linestyle = :solid, label = "33 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[1:end-219,7]),SorptionIsothermData[220:end,8]/.044, color="green", linestyle = :solid, label = "33 wt%", legend = :topleft)
plot!(sqrt.(SorptionIsothermData[1:end-279,9]),SorptionIsothermData[280:end,10]/.044, color="grey", linestyle = :solid, label = "33 wt%", legend = :topleft)

#We now fit straight lines to each of these curves. We do this individually for
#each run, so the plots are cleaner. We require each straight line to pass
#through the origin, in order for the diffusion-controlled theory to be consistent.
plot(sqrt.(SorptionIsothermData[1:end-71,1]),SorptionIsothermData[72:end,2]/.044, color="black", linestyle = :solid, xlabel = L"\sqrt{t}\mathrm{\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t}", ylims = (0,0.005), xlims = (0,140),label = "10 wt%", legend = :topleft)
plot!(LinRange(0,40,100), t -> 9.2e-5*t, color="black", linestyle = :dash, label="Diffn Contrd")

plot(sqrt.(SorptionIsothermData[1:end-319,3]),SorptionIsothermData[320:end,4]/.044, color="blue", linestyle = :solid, xlabel = L"\sqrt{t - t_0}\mathrm{\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t}", ylims = (0,0.015), xlims = (0,140),label = "29 wt%", legend = :topleft)
plot!(LinRange(0,60,100), t -> 1.55e-4*t, color="blue", linestyle = :dash, label="Diffn Contrd")

plot(sqrt.(SorptionIsothermData[1:end-379,5]),SorptionIsothermData[380:end,6]/.044, color="red", linestyle = :solid, label = "33 wt%", xlabel = L"\sqrt{t - t_0}\mathrm{\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t}", ylims = (0,0.018), xlims = (0,140), legend = :topleft)
plot!(LinRange(0,80,100), t -> 1.79e-4*t, color="red", linestyle = :dash, label="Diffn Contrd")

plot(sqrt.(SorptionIsothermData[1:end-219,7]),SorptionIsothermData[220:end,8]/.044, color="green", linestyle = :solid, label = "49 wt%", xlabel = L"\sqrt{t - t_0}\mathrm{\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t}", ylims = (0,0.025), xlims = (0,140), legend = :topleft)
plot!(LinRange(0,80,100), t -> 1.99e-4*t, color="green", linestyle = :dash, label="Diffn Contrd")

plot(sqrt.(SorptionIsothermData[1:end-279,9]),SorptionIsothermData[280:end,10]/.044, color="grey", linestyle = :solid, label = "58 wt%", xlabel = L"\sqrt{t - t_0}\mathrm{\ (s)}", ylabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t}", ylims = (0,0.025), xlims = (0,140), legend = :topleft)
plot!(LinRange(0,100,100), t -> 2.15e-4*t, color="grey", linestyle = :dash, label="Diffn Contrd")

#Same Plots with Axes Flipped
plot(SorptionIsothermData[72:end,2]/.044,sqrt.(SorptionIsothermData[1:end-71,1]), color="black", linestyle = :solid, ylabel = L"\sqrt{t-t_0}\mathrm{\ \ (s^{1/2})}", xlabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t-t_0}", xlims = (0,0.004), ylims = (0,140),label = "10 wt%", legend = :topleft)
plot!(LinRange(0,40,100), t -> t/9.2e-5, color="black", linestyle = :dash, label="Diffn Contrd")

plot(SorptionIsothermData[300:end,4]/.044,sqrt.(SorptionIsothermData[1:end-299,3]), color="blue", linestyle = :solid, ylabel = L"\sqrt{t-t_0}\mathrm{\ \ (s^{1/2})}", xlabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t - t_0}", xlims = (0,0.011), ylims = (0,140),label = "29 wt%", legend = :topleft)
plot!(LinRange(0,60,100), t -> t/1.5e-4, color="blue", linestyle = :dash, label="Diffn Contrd")

plot(SorptionIsothermData[380:end,6]/.044, sqrt.(SorptionIsothermData[1:end-379,5]),color="red", linestyle = :solid, label = "33 wt%", ylabel = L"\sqrt{t - t_0}\mathrm{\ \ (s^{1/2})}", xlabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t - t_0}", xlims = (0,0.015), ylims = (0,140), legend = :topleft)
plot!(LinRange(0,0.02,100), t -> t/1.79e-4, color="red", linestyle = :dash, label="Diffn Contrd")

plot(SorptionIsothermData[220:end,8]/.044,sqrt.(SorptionIsothermData[1:end-219,7]), color="green", linestyle = :solid, label = "49 wt%", ylabel = L"\sqrt{t - t_0}\mathrm{\ \ (s^{1/2})}", xlabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t - t_0}", xlims = (0,0.019), ylims = (0,140), legend = :topleft)
plot!(LinRange(0,80,100), t -> t/1.97e-4, color="green", linestyle = :dash, label="Diffn Contrd")

plot(SorptionIsothermData[260:end,10]/.044, sqrt.(SorptionIsothermData[1:end-259,9]), color="grey", linestyle = :solid, label = "58 wt%", ylabel = L"\sqrt{t - t_0}\mathrm{\ \ (s^{1/2})}", xlabel = L"\mathrm{mol\ CO_2}", title = L"\mathrm{CO_2\ Absorption\ vs\ } \sqrt{t - t_0}", xlims = (0,0.023), ylims = (0,140), legend = :topleft)
plot!(LinRange(0,100,100), t -> t/2.11e-4, color="grey", linestyle = :dash, label="Diffn Contrd")

#---------------------------------
#Fitting curves in a more rigorous way
#---------------------------------

#In the above section, I was fitting straight lines which passed through the
#origin to plots of $sqrt{t-t_0}$ vs $n$. However, I had to guess values of
#t_0 which allowed the line of best fit to actually pass through the orogin,
#and, frankly, it was a little dodgy. Of course, there is a much simpler and more
#rigorous way, which quite frankly I should have thought of immediately. In
#this section, we will plot

# n^2 vs t

#If the theory is correct, then we should have n^2 = C^2(t-t_0) = C^2t - C^2t_0,
#for a large part of the absorption curve (the part which is diffusion-controlled).
#From this linear region, we will be able to determine (a) C, (b) t_0, and (c)
#the range of times for which absorption was approximately diffusion controlled.
#I will do this now for each of the absorption curves above.

# 10wt% NOHMPEI
plot(SorptionIsothermData[:,1],(SorptionIsothermData[:,2]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 10wt\%\ SIP}", xlims = (0,3000), legend = :topleft,label="")

#This curve has a clear linear region. Graphically, it appears to extend from
#t = 300 to t = 1400, and so we will find a line of best fit for this region,
#and then redo the plot with this line also shown.
#We find a line of best fit using LsqFit.jl, which is called at the start of this
#file

#Define a linear model with two parameters,
@. model(x, p) = p[1]*x + p[2]
fit = curve_fit(model, SorptionIsothermData[95:407,1],(SorptionIsothermData[95:407,2]./0.044).^2,[1e-7,0.0])
fit_vals_10 = coef(fit)

#Recreate Plot with fitted values
plot(SorptionIsothermData[:,1],(SorptionIsothermData[:,2]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 10wt\%\ SIP}", xlims = (0,3000),ylims=(0,1.5e-5),label = "Experimental Data", legend = :bottomright)
plot!(LinRange(0,1550,100), t -> fit_vals_10[1]*t+fit_vals_10[2],color="black",linestyle = :dash,label="Linear Approximation")

#Extract out values of C and t0, and use these to make an n vs t plot.
C_10 = sqrt(fit_vals_10[1]); t0_10 = -fit_vals_10[2]/C_10^2
plot(SorptionIsothermData[:,1],1000*SorptionIsothermData[:,2]./0.044,color="black", linestyle = :solid,xlims=(0,4000),ylims=(0,4),legend=:bottomright,label="Experimental Data", title = L"\mathrm{CO_2\ Uptake\ in\ 10wt\%\ SIP}", xlabel = L"\mathrm{time\ (s)}", ylabel = L"\mathrm{Moles\ Absorbed\ (mmol)}")
plot!(LinRange(300,1300,15), t -> 1000*C_10*sqrt(t-t0_10),marker=:circle,color="blue",linestyle=:dash,label="Diffusion Limited Model")
plot!(LinRange(250,1500,100), t -> 1000*C_10*sqrt(t-t0_10),color="blue",linestyle=:dash,label="Diffusion Limited Model")


# 29wt% NOHMPEI
# We now do exactly the same thing for the 29wt% SIP.
plot(SorptionIsothermData[:,3],(SorptionIsothermData[:,4]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 29wt\%\ SIP}", xlims = (0,14000),label = "Experimental Data", legend = :bottomright)

#This curve has a clear linear region. Graphically, it appears to extend from
#t = 1300 to t = 3500, and so we will find a line of best fit for this region,
#and then redo the plot with this line also shown.
#We find a line of best fit using LsqFit.jl, which is called at the start of this
#file

#Fit a linear line to the data
fit = curve_fit(model, SorptionIsothermData[407:1095,3],(SorptionIsothermData[407:1095,4]./0.044).^2,[1e-7,0.0])
fit_vals_29 = coef(fit)

#Recreate Plot with fitted values
plot(SorptionIsothermData[:,3],(SorptionIsothermData[:,4]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 29wt\%\ SIP}", xlims = (0,14000),ylims=(0,0.00013),label = "Experimental Data", legend = :bottomright)
plot!(LinRange(0,5500,100), t -> fit_vals_29[1]*t+fit_vals_29[2],color="black",linestyle = :dash, label="Linear Approximation")

#Extract out values of C and t0, and use these to make an n vs t plot.
C_29 = sqrt(fit_vals_29[1]); t0_29 = -fit_vals_29[2]/C_29^2
plot(SorptionIsothermData[:,3],1000*SorptionIsothermData[:,4]./0.044,color="black", linestyle = :solid,xlims=(0,18000),ylims=(0,12),legend=:bottomright,label="Experimental Data", title = L"\mathrm{CO_2\ Uptake\ in\ 29wt\%\ SIP}", xlabel = L"\mathrm{time\ (s)}", ylabel = L"\mathrm{Moles\ Absorbed\ (mmol)}",bottom_margin = 5mm)
plot!(LinRange(1100,3500,15), t -> 1000*C_29*sqrt(t-t0_29),marker=:circle,color="blue",linestyle=:dash,label="Diffusion Limited Model")
plot!(LinRange(950,5000,100), t -> 1000*C_29*sqrt(t-t0_29),color="blue",linestyle=:dash,label="Diffusion Limited Model")

# 33wt% NOHMPEI
# We now do exactly the same thing for the 33wt% SIP.
plot(SorptionIsothermData[:,5],(SorptionIsothermData[:,6]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t \ \mathrm{for\ 33wt\%\ SIP}" ,xlims = (0,14000),ylims=(0,0.00025),label = "Experimental Data", legend = :bottomright,bottom_margin = 5mm)

#This curve has a clear linear region. Graphically, it appears to extend from
#t = 1500 to t = 5500, and so we will find a line of best fit for this region,
#and then redo the plot with this line also shown.
#We find a line of best fit using LsqFit.jl, which is called at the start of this
#file

#Fit a linear line to the data
fit = curve_fit(model, SorptionIsothermData[470:1720,5],(SorptionIsothermData[470:1720,6]./0.044).^2,[1e-7,0.0])
fit_vals_33 = coef(fit)

#Recreate Plot with fitted values
plot(SorptionIsothermData[:,5],(SorptionIsothermData[:,6]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 33wt\%\ SIP}", xlims = (0,14000),ylims=(0,0.00025),label = "Experimental Data", legend = :bottomright,bottom_margin = 5mm)
plot!(LinRange(0,7000,100), t -> fit_vals_33[1]*t+fit_vals_33[2],color="black",linestyle = :dash, label = "Linear Approximation")

#Extract out values of C and t0, and use these to make an n vs t plot.
C_33 = sqrt(fit_vals_33[1]); t0_33 = -fit_vals_33[2]/C_33^2
plot(SorptionIsothermData[:,5],1000*SorptionIsothermData[:,6]./0.044,color="black", linestyle = :solid,xlims=(0,18000),ylims=(0,15),legend=:bottomright,label="Experimental Data", title = L"\mathrm{CO_2\ Uptake\ in\ 33wt\%\ SIP}", xlabel = L"\mathrm{time\ (s)}", ylabel = L"\mathrm{Moles\ Absorbed\ (mmol)}",bottom_margin = 5mm)
plot!(LinRange(1500,5500,15), t -> 1000*C_33*sqrt(t-t0_33),color="blue",marker=:circle,linestyle=:dash,label="Diffusion Limited Model")
plot!(LinRange(1200,7000,100), t -> 1000*C_33*sqrt(t-t0_33),color="blue",linestyle=:dash,label = "Diffusion Limited Model")

# 49wt% NOHMPEI
# We now do exactly the same thing for the 49wt% SIP.
plot(SorptionIsothermData[:,7],(SorptionIsothermData[:,8]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 49wt\%\ SIP}", xlims = (0,18000),ylims=(0,0.00045),label = "", legend = :topleft,bottom_margin = 5mm)

#The linear region is less obvious here. Graphically, it roughly extends from
#t = 1500 to t = 5000, and so we will find a line of best fit for this region,
#and then redo the plot with this line also shown.
#We find a line of best fit using LsqFit.jl, which is called at the start of this
#file

#Fit a linear line to the data
fit = curve_fit(model, SorptionIsothermData[470:1566,7],(SorptionIsothermData[470:1566,8]./0.044).^2,[1e-7,0.0])
fit_vals_49 = coef(fit)

#Recreate Plot with fitted values
plot(SorptionIsothermData[:,7],(SorptionIsothermData[:,8]/.044).^2, color="black", linestyle = :solid, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 49wt\%\ SIP}", xlabel = L"\mathrm{time\ (s)}", xlims = (0,18000),ylims=(0,0.0004),label = "Experimental Data", legend = :bottomright, bottom_margin = 5mm)
plot!(LinRange(0,9000,100), t -> fit_vals_49[1]*t+fit_vals_49[2],color="black",linestyle = :dash,label = "Linear Approximation")

#Extract out values of C and t0, and use these to make an n vs t plot.
C_49 = sqrt(fit_vals_49[1]); t0_49 = -fit_vals_49[2]/C_49^2
plot(SorptionIsothermData[:,7],1000*SorptionIsothermData[:,8]./0.044,color="black", linestyle = :solid,xlims=(0,18000),ylims=(0,20),legend=:bottomright,label="Experimental Data", title = L"\mathrm{CO_2\ Uptake\ in\ 49\%\ SIP}", xlabel = L"\mathrm{time\ (s)}", ylabel = L"\mathrm{Moles\ Absorbed\ (mmol)}",bottom_margin = 5mm)
plot!(LinRange(1500,5000,15), t -> 1000*C_49*sqrt(t-t0_49),color="blue",marker=:circle,linestyle=:dash,label = "Diffusion Limited Model")
plot!(LinRange(600,7500,15), t -> 1000*C_49*sqrt(t-t0_49),color="blue",linestyle=:dash,label = "Diffusion Limited Model")


# 58wt% NOHMPEI
# We now do exactly the same thing for the 58wt% SIP.
plot(SorptionIsothermData[:,9],(SorptionIsothermData[:,10]/.044).^2, color="black", linestyle = :solid,  title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 58wt\%\ SIP}", xlims = (0,18000),ylims=(0,0.00068),label = "", legend = :topleft)

#The linear region is less obvious here. Graphically, it roughly extends from
#t = 1500 to t = 5000, and so we will find a line of best fit for this region,
#and then redo the plot with this line also shown.
#We find a line of best fit using LsqFit.jl, which is called at the start of this
#file

#Fit a linear line to the data
fit = curve_fit(model, SorptionIsothermData[470:1566,9],(SorptionIsothermData[470:1566,10]./0.044).^2,[1e-7,0.0])
fit_vals_58 = coef(fit)

#Recreate Plot with fitted values
plot(SorptionIsothermData[:,9],(SorptionIsothermData[:,10]/.044).^2, color="black", linestyle = :solid,xlabel=L"\mathrm{time\ (s)}",bottom_margin=5mm, title = L"n^2\ \mathrm{vs}\ t\ \mathrm{for\ 58wt\%\ SIP}", xlims = (0,18000),ylims=(0,0.0006),label = "Experimental Data", legend = :bottomright)
plot!(LinRange(0,9000,100), t -> fit_vals_58[1]*t+fit_vals_58[2],color="black",linestyle = :dash,label="Linear Approximation")

#Extract out values of C and t0, and use these to make an n vs t plot.
C_58 = sqrt(fit_vals_58[1]); t0_58 = -fit_vals_58[2]/C_58^2
plot(SorptionIsothermData[:,9],1000*SorptionIsothermData[:,10]./0.044,color="black", linestyle = :solid,xlims=(0,18000),ylims=(0,25),legend=:bottomright,label="Experimental Data", title = L"\mathrm{CO_2\ Uptake\ in\ 58wt\%\ SIP}", xlabel = L"\mathrm{time\ (s)}", ylabel = L"\mathrm{Moles\ Absorbed\ (mmol)}",bottom_margin = 5mm)
plot!(LinRange(1500,5000,15), t -> 1000*C_58*sqrt(t-t0_58),color="blue",marker=:circle,linestyle=:dash,label = "Diffusion Limited Model")
plot!(LinRange(860,9000,100), t -> 1000*C_58*sqrt(t-t0_58),color="blue",linestyle=:dash, label = "Diffusion Limited Model")


#---------------------------------------------------------#
#Calculation of Predicted SIP Thicknesses from Uptake Data
#---------------------------------------------------------#

#From the uptake curves, we have identified regions in which absorption is likely
#to be diffusion controlled. We can now back-calculate from the absorption data
#to determine the predicted thickness of the SIP layer. If the theory is correct,
#this should correspond to a value between 0.5mm (double-sided abs.) and 1mm
#(single-sided absorption). We will require:

# - CO2 Partial Pressure, pCO2 = 1 bar = 100,000 Pa
# - Permeability of CO2 in PDMS saturated with N2 = 1650 barrer = 5.58e-13 mol/m.Pa.s (Scholes et al. J. Mem. Sci. 2010)
# - Solubility within PDMS = 5.94e-4 mol/Pa.m3
# - Capacity within NPEI = 5.5 mol/kg
# - We set the density of NPEI at 1250 kg/m3, however this cancels during the calculation, so the exact value is irrelevant.

SPDMS = 5.94e-4                                             #mol/m3.Pa
ρPDMS = 965                                                 #kg/m3
PCO2 = 1e5                                                  #Pa
pCO2 = 1e5                                                      #Pressure of CO2 Pa
PPDMS = 5.58e-13                                                #CO2 Permeability of PDMS mol/m.Pa.s
PPDMS = 1.1e-12
PNPEI = 0                                                       #CO2 Permeability of NPEI mol/m.Pa.s
CapNPEI = 5.5                                                   #CO2 Capacity of NPEI mol/kg
ρNPEI = 1250                                                    #Density of NPEI kg/m3

#Beging calculations of theoretical lengths:
εMass =  [1-MassFractions[i] for i in 1:5]                                     #Mass Fraction of PDMS
ε = εMass./ρPDMS ./ (εMass./ρPDMS .+ (1 .-εMass)./ρNPEI)
N = ε.*SPDMS.*pCO2 .+ (1 .-ε).*ρNPEI.*CapNPEI

#Now I need to calculate C1 values from the C values which were calculated numerically.
#Note that C and C1 are not the same - see Eq. 4.34 in thesis for the definition
#of C1.

#We now calculate the total theoretical capacity of each SIP material.
N_max_SIPs = [SIPMasses[i]*(1-εMass[i])*CapNPEI + SIPMasses[i]*εMass[i]/ρPDMS*SPDMS*pCO2 for i in 1:5]
C_data = [C_10,C_29,C_33,C_49,C_58]
C1_data = [(C_data[i]/N_max_SIPs[i])^2 for i in 1:5]
L = sqrt.(2 .*pCO2./(N.*C1_data).*((1 .-ε).*PNPEI .+ ε.*PPDMS))
println(L)
