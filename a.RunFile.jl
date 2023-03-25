#RUN FILE

# Calls the Packages used for the optimization problem
using JuMP
using Printf
using CPLEX
using MathOptInterface
using JLD
using TimerOutputs
using Distributions
using DataFrames
using XLSX
using Parameters
using Dates
using CSV
using Plots
using Combinatorics
using Rainflow
using Base

import Base.show

# Calls the other Julia files
include("Structures.jl")
include("SetInputParameters.jl")
include("dynamicProgramming.jl")
include("Saving in xlsx.jl")
include("Plots.jl")

date = string(today())

# PREPARE INPUT DATA
to = TimerOutput()

@timeit to "Set input data" begin

  #Set run case - indirizzi delle cartelle di input ed output
  case = set_runCase()
  @unpack (DataPath,InputPath,ResultPath,CaseName) = case;

  # Set run mode (how and what to run) and Input parameters
  runMode = read_runMode_file()
  InputParameters = set_parameters(runMode, case)
  @unpack (NStages, NStates, NHoursStage, Big)= InputParameters;

  # Set solver parameters (Cplex etc)
  SolverParameters = set_solverParameters()

  # Read power prices from a file [€/MWh]
  priceYear = "2022"
  Power_prices = read_csv(" Prezzi mezz'ora 2022.csv",case.DataPath)                       # valori alla mezz'ora
  #Battery_prices = read_csv("Cost_battery.csv",case.DataPath)                             #  cost for battery replacement for each half hour
  PV_production = read_csv("PV_production.csv",case.DataPath)                              # potenza da PV (0.75MW)         

  # Upload battery's characteristics
  Battery = set_battery_system(runMode, case)
  @unpack (power_Capacity, energy_Capacity, Eff_charge, Eff_discharge, costBattery, DoD, NCycles) = Battery; 

  # DEFINE STATE VARIABLES - STATE OF CHARGES SOC [MWh]
  state_variables = define_state_variables(InputParameters, Battery)

  # Where and how to save the results
  FinalResPath= set_run_name(case, ResultPath, InputParameters)

end

#save input data
@timeit to "Save input" begin
    save(joinpath(FinalResPath,"CaseDetails.jld"), "case" ,case)
    save(joinpath(FinalResPath,"SolverParameters.jld"), "SolverParameters" ,SolverParameters)
    save(joinpath(FinalResPath,"InputParameters.jld"), "InputParameters" ,InputParameters)
    save(joinpath(FinalResPath,"BatteryCharacteristics.jld"), "BatteryCharacteristics" ,Battery)
    save(joinpath(FinalResPath,"PowerPrices.jld"),"PowerPrices",Power_prices)
    save(joinpath(FinalResPath,"PVproduction.jld"),"PVprod",PV_production)
end

configurations=[]
push!(configurations, " Stand-alone system BESS")
push!(configurations, " BESS with PV")
push!(configurations, " BESS with PV and maximum Grid capacity")
push!(configurations, " BESS with PV - only discharge")

# DYNAMIC PROGRAMMING
if runMode.dynamicProgramming
    println("Solving Dynamic Pogramming")
    ResultsDP = DP(InputParameters, Battery, state_variables, runMode, Power_prices, PV_production)   #configurations
    save(joinpath(FinalResPath, "dp_Results.jld"), "dp_Results", ResultsDP) 
    else
    println("Solved without dynamic programming.")
end

# SAVE OTIMAL-PATH DATA IN EXCEL FILES
if runMode.excel_savings
  cartella = "C:\\GitSource-Batteries\\Batteries-greedy-aproach\\Results"
  cd(cartella)
  data_saving(InputParameters,ResultsDP,Power_prices)
  println("Results saved")
else
  println("Solved without saving results in xlsx format.")
end

# SAVE PLOTS IN THE CORRESPONDING FOLDER
if runMode.plot_savings
  cartella = "C:\\GitSource-Batteries\\Batteries-greedy-aproach\\Plots"
  cd(cartella)
  plotPath(InputParameters,ResultsDP,state_variables,priceYear)
end

print(to)




