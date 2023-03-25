function DP(                                                                                                   # Ora conosco per ogni settimana i valori di inflow e prezzo (calcolati con il modello di Markov) - risolvo il problema come "DETERMINISTICO"
  InputParameters::InputParam,
  Battery::BatteryParam,
  state_variables::states,
  runMode::runModeParam,
  Power_prices,
  PV_production,
  )

  @unpack (NStages, NStates, NHoursStage, Big) = InputParameters
  @unpack (grid_Capacity, power_Capacity, energy_Capacity, Eff_charge, Eff_discharge, costBattery, DoD, NCycles) = Battery      # MAXCharge, MAXDischarge,
 # @unpack (CPX_PARAM_SCRIND, CPX_PARAM_PREIND, CPXPARAM_MIP_Tolerances_MIPGap, CPX_PARAM_TILIM, CPX_PARAM_THREADS) = SolverParameters
  @unpack (productionPV, battery_replacement) = runMode
  @unpack (seg) = state_variables

  power_price = Power_prices[1:NStages]
  pv_prod = zeros(NStages)
  costBat = 0

  mode = "Optimal trajectory"
  batCost=" "
  if battery_replacement
    batCost="_with Battery replacement"
    costBat = costBattery
  else
    batCost="_without Battery replacement"
  end

  Pvprod=""
  if productionPV
    Pvprod="_with PV production"
    pv_prod= PV_production[1:NStages]
  else
    Pvprod="_without PV production"
  end

  modeFinal = mode*batCost*Pvprod
  println(modeFinal)

  optimalValueStates = zeros(NStages+1,NStates)                                # Final Optimal value for each State of the t-stage -> considers the max among all combinations ex: ValueStates[23,5] = 124€ -> if we are in day 23 at stage 5, the value I would have is of 124€
  optimalValueStates[end,:] = seg * Power_prices[NStages+1]
  optimalfromState = zeros(NStages,NStates)                                    # Indicates the optimal state of the next stage from which we are coming from ex: fromState[23,5] =2 -> if we are at day 24 in state 5 (0% of energy), we are comiing from state 2 in day 24
  val = zeros(NStages,NStates,NStates)                                         # Per ogni stato del sistema, calcolo tutte le transizioni e poi ne prendo la massima                                                              # ex: Val[35,1,4] = indicates the value of being in state 1 in day 35, id coming from state 4 in stage 36

  # VECTORS FOR EVERY POSSIBLE COMBINATION
  charge = zeros(NStages,NStates,NStates)
  charge_grid = zeros(NStages,NStates,NStates) 
  discharge = zeros(NStages,NStates,NStates)
  degradation = zeros(NStages,NStates,NStates)
  optimalPath = []
  gain = zeros(NStages,NStates,NStates)
  replacementCost = zeros(NStages,NStates,NStates)

  @timeit to "Solve dynamic programming" begin

    for t = NStages:-1:1                                                        # Calcolo per ogni mezz'ora partendo dall'ultimo
      println("t:", t)
      
      soc_start = 0
      soc_final = 0
      penalty = 0

      for iState=1:NStates                                                      # Considero gg=365gg*24h*2 tutti i possibili stati

        soc_start = seg[iState]

          for jState=1:NStates                                                  # Considero tutti gli stati allo stagesuccessivo

            #CALCULATES THE CHARGE/DISCHARGE FROM ONE STAGE TO ANOTHER CONSIDERING ALL POSSIBLE STATE TRANSITIONS PER EACH STAGE

            soc_final = seg[jState]

            if soc_final > soc_start
                charge[t,iState,jState] = abs((soc_final-soc_start)/(NHoursStage*Eff_charge))         #calculates how much power is needed to charge the battery from one stato to another
                discharge[t,iState,jState]=0
                degradation[t,iState,jState]=abs(1/NCycles[iState]*0.5-1/NCycles[jState]*0.5)

                # INFEASIBILITIES FOR CHARGE>2MW - high cost associated
                if charge[t,iState,jState]>power_Capacity
                  penalty = Big
                else
                  penalty = 0
                end
                
                # CONSIDERS CONTRIBUTION FROM PV
                if charge[t,iState,jState]<pv_prod[t]
                  charge_grid[t,iState,jState] = 0                  
                else
                  charge_grid[t,iState,jState] = charge[t,iState,jState]-pv_prod[t]
                end
                
            elseif soc_final < soc_start
                charge[t,iState,jState] = 0
                charge_grid[t,iState,jState] = 0
                discharge[t,iState,jState] = abs((soc_final-soc_start)*Eff_discharge/NHoursStage)
                degradation[t,iState,jState] = abs(1/NCycles[iState]*0.5-1/NCycles[jState]*0.5)              #*10E5
                penalty = 0
            else
                charge[t,iState,jState] = 0
                charge_grid[t,iState,jState] = 0
                discharge[t,iState,jState] = 0
                degradation[t,iState,jState] = 0
                penalty = 0
            end

           
            val[t,iState,jState] = power_price[t]*NHoursStage*(discharge[t,iState,jState]-charge_grid[t,iState,jState]) - degradation[t,iState,jState]*costBat*energy_Capacity - penalty + optimalValueStates[t+1,jState]      #/10E5
            gain[t,iState,jState] = power_price[t]*NHoursStage*(discharge[t,iState,jState]-charge_grid[t,iState,jState])
            replacementCost[t,iState,jState] = degradation[t,iState,jState]*costBat*energy_Capacity                                                                 #/10E5

          end # end jStates=1:5

        optimalValueStates[t,iState] = findmax(val[t,iState,:])[1]             # Trovo il massimo del Valore funzione obiettivo : transizioni + valore stato precedente 
        optimalfromState[t,iState] = findmax(val[t,iState,:])[2]               # Mi dice da quale stato al giorno precedente (o futuro) arrivo

      end

    end   # end Stages=1:365

    # RACCOLGO I RISULTATI DEL PERCORSO MIGLIORE

    #startingValue = findmax(optimalValueStates[1,:])[1]                                 # Indica il valore ottimo (massimo Val) allo Stage = 1 
    startingFrom = findmax(optimalValueStates[1,:])[2]                                  # Mi inidica in quale stato per NStage=1 ho il massimo valore
    netOverallGain = 0
    overallCost = 0

    for t=1:NStages   
      
      comingFrom = Int(optimalfromState[t,startingFrom])                                       # Inidca da quale stato presedente sono arrivato
    
      optValue = findmax(optimalValueStates[t,startingFrom])[1]
      final_operation=charge[t,startingFrom,comingFrom]-discharge[t,startingFrom,comingFrom]
      tot_gain = gain[t,startingFrom,comingFrom] - replacementCost[t,startingFrom,comingFrom]
      batteryCost = replacementCost[t,startingFrom,comingFrom]  

      overallCost = overallCost+batteryCost
      netOverallGain = netOverallGain+tot_gain

      push!(optimalPath,saveOptimalValues(t,optValue,startingFrom,comingFrom,seg,final_operation,tot_gain,batteryCost))
      startingFrom=comingFrom

    end
    
    error = []
    for t=1:NStages 
      if optimalPath[t].action> power_Capacity
        push!(error,optimalPath[t])
      end
    end

  end

  return Results_dp(
    charge,
    charge_grid,
    discharge,
    degradation,
    gain,
    pv_prod,
    replacementCost,
    power_price,
    val,
    optimalValueStates,
    optimalfromState,
    optimalPath,
    overallCost,
    netOverallGain,
    error,
    modeFinal,
   )

end

show(io::IO,x::optimalStage) = print(io, "Stage:",x.stage," -> optimal Value = ",x.optimalValue," ,current State = ",x.currentState," ,next State = ",x.nextState," ,current SOC = ",x.currentSOC," ,next SOC = ",x.nextSOC," ,action = ",x.action, " , net gain =",x.gain, "battery Cost =",x.batteryCost)

function saveOptimalValues(stage::Int64,optimalValue::Float64, curSt::Int64, nextSt::Int64, seg::Any, decision::Float64, gain::Float64 , batteryCost::Float64)
  optimalStage(stage,optimalValue,curSt,nextSt,seg[curSt],seg[nextSt],decision,gain,batteryCost)
end