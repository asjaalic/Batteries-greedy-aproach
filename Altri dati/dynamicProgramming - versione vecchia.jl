function DP(                                                                                                   # Ora conosco per ogni settimana i valori di inflow e prezzo (calcolati con il modello di Markov) - risolvo il problema come "DETERMINISTICO"
  InputParameters::InputParam,
  Battery::BatteryParam,
  state_variables::states,
  runMode::runModeParam,
  Power_prices,
  PV_production,
  )

  @unpack (NStages, NStates, NHoursStage, Big) = InputParameters
  @unpack (grid_Capacity, power_Capacity, energy_Capacity, Eff_charge, Eff_discharge,costBattery, DoD, NCycles) = Battery      # MAXCharge, MAXDischarge,
  @unpack (maximumGrid, onlyExport, productionPV, battery_replacement) = runMode
  @unpack (seg) = state_variables

  power_price = Power_prices[1:NStages]
  pv_prod = zeros(NStages)
  costBat = 0
  gridMax = grid_Capacity

  mode = "Optimal trajectory"

  batCost=" "
  if battery_replacement                                                        # if true -> evaluating with cell -replacement , give the cost
    batCost="_with Battery replacement"
    costBat = costBattery
  else
    batCost="_without Battery replacement"
  end

  Pvprod=""
  if productionPV                                                               # if true -> consider the production from PV pannels
    Pvprod="_with PV production"
    pv_prod= PV_production[1:NStages]
  else
    Pvprod="_without PV production"
  end

  modeFinal = mode*batCost*Pvprod                                               # full name of the problem solving
  println(modeFinal)
  
  optimalValueStates = zeros(NStages+1,NStates)                                 # Final Optimal value for each State of the t-stage -> considers the max among all combinations ex: ValueStates[23,5] = 124€ -> if we are in day 23 at stage 5, the value I would have is of 124€
  optimalValueStates[end,:] = seg * Power_prices[NStages+1]                     # Initialize the Values of NStages+1 (starting point DP)
  optimalfromState = zeros(NStages,NStates)                                     # Indicates the optimal state from which we are coming from ex: fromState[23,5] =2 -> if we are at day 23 in state 5 (0% of energy), we are comiing from state 2 in day 24
  val = zeros(NStages,NStates,NStates)                                          # Per ogni stato del sistema, calcolo tutte le transizioni e poi ne prendo la massima                                                              # ex: Val[35,1,4] = indicates the value of being in state 1 in day 35, id coming from state 4 in stage 36
 

  # VECTORS FOR EVERY POSSIBLE COMBINATION
  charge = zeros(NStages,NStates,NStates)                                       # power needed to charge the battery
  discharge = zeros(NStages,NStates,NStates)                                    # power discharged by the battery
  discharge_grid = zeros(NStages,NStates,NStates)                               # effective power needed from the grid
  charge_grid =zeros(NStages,NStates,NStates)                                   # effective power sold to the grid

  degradation = zeros(NStages,NStates,NStates)                                  # accounts for the %of battery degradated beacuse of the use
  optimalPath = []

  gain = zeros(NStages,NStates,NStates)
  replacementCost = zeros(NStages,NStates,NStates)

  @timeit to "Solve dynamic programming" begin

    for t = NStages:-1:1                                                        # Calcolo per ogni mezz'ora partendo dall'ultimo
      println("t:", t)
      
      soc_start = 0
      soc_final = 0
      penaltySOC = 0
      penaltyExport = 0

      for iState=1:NStates                                                      # Considero gg=365gg*24h*2 tutti i possibili stati

        soc_start = seg[iState]

          for jState=1:NStates                                                  # Considero tutti gli stati allo stagesuccessivo

            #CALCULATES THE CHARGE/DISCHARGE FROM ONE STAGE TO ANOTHER CONSIDERING ALL POSSIBLE STATE TRANSITIONS PER EACH STAGE

            soc_final = seg[jState]
            #charge_grid[t,iState,jState]= charge[t,iState,jState]-pv_prod[t]
            #discharge_grid[t,iState,jState] = discharge[t,iState,jState]+pv_prod[t]

            if soc_final > soc_start          #CHARGING PHASE
                
                charge[t,iState,jState] = abs((soc_final-soc_start)/(NHoursStage*Eff_charge))         #calculates how much power is needed to charge the battery from one stato to another
                discharge[t,iState,jState] = 0
                degradation[t,iState,jState]=abs(1/NCycles[iState]*0.5-1/NCycles[jState]*0.5)

                charge_grid[t,iState,jState] = charge[t,iState,jState]-pv_prod[t]
                
                # INFEASIBILITIES FOR CHARGE>2MW
                if charge[t,iState,jState]>power_Capacity
                  penaltySOC = Big
                end
                
                if onlyExport                 # se posso solo caricare la batteria solo dal PV, non datta rete  
                  if charge_grid > 0          # se devo prendere dalla rete, metto una penale
                    penaltyExport = Big
                    discharge_grid[t,iState,jState] = 0
                  else                        # se ho più PV di quanto ne abbia bisogno per caricare la batteria
                    charge_grid[t,iState,jState] = 0
                    discharge_grid[t,iState,jState] = pv_prod[t]-charge[t,iState,jState]
                  end
                else                          # false -> se posso caricare la batteria anche dalla rete oltre che da PV   
                  if charge[t,iState,jState]<=pv_prod[t]
                    charge_grid[t,iState,jState] = 0 
                    discharge_grid[t,iState,jState] = pv_prod[t]-charge[t,iState,jState]                 
                  else                        # otherwise what is missing, have to buy from the grid
                    charge_grid[t,iState,jState] = charge[t,iState,jState]-pv_prod[t]
                    discharge_grid[t,iState,jState] = 0
                  end

                end
                   
            elseif soc_final < soc_start       #DISCHARGING PHASE
              
              charge[t,iState,jState] = 0
              discharge[t,iState,jState] = abs((soc_final-soc_start)*Eff_discharge/NHoursStage)       #discharge from battery
              degradation[t,iState,jState] = abs(1/NCycles[iState]*0.5-1/NCycles[jState]*0.5)         

              charge_grid[t,iState,jState] = 0
              discharge_grid[t,iState,jState] = discharge[t,iState,jState]+pv_prod[t]

              # se vi è un limite sulla capacità di rete e viene superato ,scarico gridMax - altrimenti scarico tutto
              if maximumGrid && discharge_grid[t,iState,jState]>=gridMax
                discharge_grid[t,iState,jState] = gridMax
              end
                  
            else  #se non c'è nè carica nè scarica, posso solo vendere l'energia dei PV alla rete
                
              charge[t,iState,jState] = 0
              discharge[t,iState,jState] = 0
              degradation[t,iState,jState] = 0

              charge_grid[t,iState,jState] = 0
              discharge_grid[t,iState,jState]=discharge[t,iState,jState]+pv_prod[t]

              if maximumGrid && discharge_grid[t,iState,jState]>= gridMax #se vi è un limite su potenza massima nella rete e la potenza da PV è maggiore della capacità di rete)
                  discharge_grid[t,iState,jState] = gridMax
              end
      
            end

            val[t,iState,jState] = power_price[t]*NHoursStage*(discharge_grid[t,iState,jState]-charge_grid[t,iState,jState]) - degradation[t,iState,jState]*costBattery*energy_Capacity - penaltySOC - penaltyExport + optimalValueStates[t+1,jState]      #/10E5
            gain[t,iState,jState] = power_price[t]*NHoursStage*(discharge_grid[t,iState,jState]-charge_grid[t,iState,jState])
            replacementCost[t,iState,jState] = degradation[t,iState,jState]*costBattery*energy_Capacity                                                                 #/10E5

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
      charge_bat = charge[t,startingFrom,comingFrom]
      dis_bat = discharge[t,startingFrom,comingFrom]
      char_grid = charge_grid[t,startingFrom,comingFrom]
      dis_grid = discharge_grid[t,startingFrom,comingFrom]
      pv = pv_prod[t]

      tot_gain = gain[t,startingFrom,comingFrom] - replacementCost[t,startingFrom,comingFrom]
      batCost = replacementCost[t,startingFrom,comingFrom]  

      overallCost = overallCost+batCost
      netOverallGain = netOverallGain+tot_gain

      push!(optimalPath,saveOptimalValues(t,optValue,startingFrom,comingFrom,seg,charge_bat,char_grid,dis_bat,dis_grid,pv,tot_gain,batCost))
      startingFrom=comingFrom

    end   
  end

  return Results_dp(
    charge,
    discharge,
    charge_grid,
    discharge_grid,
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


show(io::IO,x::optimalStage) = print(io, "Stage:",x.stage," -> opVal = ",x.optimalValue," ,curSOC = ",x.currentSOC," ,nextSOC = ",x.nextSOC," ,charge = ",x.charge, ",fromGrid = ",x.charge_grid, "discharge = ",x.discharge , "toGrid = ",x.discharge_grid,"PV =", x.PV , "net gain =",x.gain, "battery Cost =",x.batteryCost)

function saveOptimalValues(stage::Int64,optimalValue::Float64, curSt::Int64, nextSt::Int64, seg::Any, charge::Float64,discharge::Float64,fromGrid::Float64, toGrid::Float64, PV::Float64, gain::Float64 , batteryCost::Float64)
  optimalStage(stage,optimalValue,seg[curSt],seg[nextSt],charge,fromGrid,discharge,toGrid,PV,gain,batteryCost)
end



#=
   tot = discharge[t,iState,jState]+pv_prod[t]

              if maximumGrid  #(true - se vi è un limite su potenza massima nella rete e il totale lo supera)
                if tot>=gridMax 
                  discharge_grid[t,iState,jState] = gridMax
                else      # se c'è un limite ma non viene superato, vendo quello che ho
                  discharge_grid[t,iState,jState] = discharge[t,iState,jState]+pv_prod[t]
                end
              else #(false - non vi è un limite sulla potenza massima nella rete -> vendo tutto quello che ho)
                discharge_grid[t,iState,jState] = discharge[t,iState,jState]+pv_prod[t]
              end
              
              
              =#