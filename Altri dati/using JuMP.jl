using JuMP
using CPLEX
NSteps=48

M = Model(
      with_optimizer(
        CPLEX.Optimizer,
      ),
    )

@variable(M,xc,Bin, base_name = "Charging")
@variable(M,xd,Bin, base_name = "Discharging")

@variable(M, 0 <= soc[iStep=1:NSteps] <= 10, base_name = "State of Charge")                # MWh   power_capacity
@variable(M, 0 <= charge[iStep=1:NSteps] <= 10, base_name = "Charge")                       # MAXCharge  - variabile decisionale : ricarico la batteria MW
@variable(M, 0 <= discharge[iStep=1:NSteps] <= 10, base_name = "Discharge")                 # MAXDischarge - variabile decisionale: scarico la batteria
@variable(M, value_next <= 10000, base_name = "Value next stage")

@objective(M,MathOptInterface.MAX_SENSE, sum(1*discharge[iStep]-1*charge[iStep] for iStep=1:NSteps)+value_next)
         
# UPDATE CONSTRAINTS

@constraint(M,Final_state_variable[iStep=NSteps], soc[iStep] + (Eff_charge*charge[iStep]-Eff_discharge*discharge[iStep])*NHoursStep == 0)
@constraint(M,State_variable[iStep=1:(NSteps-1)], soc[iStep] + (Eff_charge*charge[iStep]-Eff_discharge*discharge[iStep])*NHoursStep == soc[iStep+1])
@constraint(M,endingE[iStep=1], soc[iStep] == energy_Capacity)
@constraint(M,Value_next_stage, value_next == energy_Capacity)

@constraint(M, Binary_variables, xc+xd==1)
@constraint(M,Charging[iStep=1:NSteps],charge[iStep] <= 10*xc)
@constraint(M,Discharging[iStep=1:NSteps], discharge[iStep]<=10*xd)