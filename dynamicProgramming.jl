function DP(                                                                                                   # Ora conosco per ogni settimana i valori di inflow e prezzo (calcolati con il modello di Markov) - risolvo il problema come "DETERMINISTICO"
  InputParameters::InputParam,
  Battery::BatteryParam,
  state_variables::states,
  runMode::runModeParam,
  Power_prices,
  PV_production,
  #configurations,
  )

  @unpack (NStages, NStates, NHoursStage, Big) = InputParameters
  @unpack (grid_Capacity, power_Capacity, energy_Capacity, Eff_charge, Eff_discharge,costBattery, DoD, NCycles) = Battery      # MAXCharge, MAXDischarge,
  @unpack (maximumGrid, onlyExport, productionPV, battery_replacement) = runMode
  @unpack (seg) = state_variables

  power_price = Power_prices[1:NStages]       # 1:NStages
  pv_prod = zeros(NStages)
  costBat = 0
  gridMax = grid_Capacity #0

  bat_Cost=" "
  if battery_replacement                                                        # if true -> evaluating with cell -replacement , give the cost
    bat_Cost="_WITHBatRep"
    costBat = costBattery
  else
    bat_Cost="_WITHOUTBatRep"
  end

  #ATT - maximumGrid e only Export cannot be true at same time!!

  modePV=""
  if productionPV
    modePV= "_withPV"
    pv_prod = PV_production[1:NStages]
  else
    modePV= "_withoutPV"
  end

  modeGrid=""
  if maximumGrid
    modeGrid= "_withMaxGrid"
  else
    modeGrid= "_withoutMaxGrid"
  end

  modeExport= ""
  if onlyExport
    modeExport= "_onlyExport"
  else
    modeExport= "_Imp&Exp"
  end

  #=
  if productionPV && maximumGrid                                                               # if true -> consider the production from PV pannels
    mode=configurations[3]
    pv_prod= PV_production[1:NStages]
    #gridMax = grid_Capacity
  elseif productionPV && !maximumGrid
    mode=configurations[2]
    pv_prod = PV_production[1:NStages]
  elseif productionPV && onlyExport
    mode = configurations[4]
  elseif !productionPV
    mode=configurations[1]
    #pv_prod=zeros(NStages)
  end =#

  modeFinal = bat_Cost*modePV*modeGrid*modeExport                                                      # full name of the problem solving
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
      println("STAGE:", t)
      
      soc_start = 0
      soc_final = 0

      for iState=1:NStates                                                      # Considero gg=365gg*24h*2 tutti i possibili stati

        soc_start = seg[iState]

          for jState=1:NStates                                                  # Considero tutti gli stati allo stagesuccessivo

            #CALCULATES THE CHARGE/DISCHARGE FROM ONE STAGE TO ANOTHER CONSIDERING ALL POSSIBLE STATE TRANSITIONS PER EACH STAGE

            soc_final = seg[jState]
            penaltySOC = 0
            penaltyExport = 0

            if soc_final > soc_start          #CHARGING PHASE
                
                charge[t,iState,jState] = abs((soc_final-soc_start)/(NHoursStage*Eff_charge))         #calculates how much power is needed to charge the battery from one stato to another
                discharge[t,iState,jState] = 0
                degradation[t,iState,jState]=abs(1/NCycles[iState]*0.5-1/NCycles[jState]*0.5)

                charge_grid[t,iState,jState] = charge[t,iState,jState]-pv_prod[t]
                discharge_grid[t,iState,jState] = 0        #discharge[t,iState,jState] + pv_prod[t]
                
                # INFEASIBILITIES FOR CHARGE>POwER_CAPACITY MW - if true add penalty, otherwise leave 0
                if charge[t,iState,jState]>power_Capacity
                  penaltySOC = Big
                else
                  penaltySOC = 0
                end
                
                if onlyExport                                                                         # se posso caricare la batteria solo dal PV, non daa rete  
                  if charge_grid[t,iState,jState] > 0                                                                  # se devo prendere dalla rete, metto una penale
                    penaltyExport = 10E7
                    #discharge_grid[t,iState,jState] = 0
                    #charge_grid[t,iState,jState] = 0
                  else
                    penaltyExport = 0                                                                                # se ho più PV di quanto ne abbia bisogno per caricare la batteria
                    charge_grid[t,iState,jState] = 0
                    discharge_grid[t,iState,jState] = pv_prod[t]-charge[t,iState,jState]
                  end
                else
                  penaltyExport = 0                                                                                 # false -> se posso caricare la batteria anche dalla rete oltre che da PV   
                  if charge_grid[t,iState,jState] <= 0                                                # se ho più PV di quanto ne abbia bisogno per caricare la batteria, vendo il restante in rete
                    charge_grid[t,iState,jState] = 0 
                    discharge_grid[t,iState,jState] = pv_prod[t]-charge[t,iState,jState]                 
                  else                                                                                # otherwise what is missing, have to buy from the grid
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
                #penaltyExport = 10 ??
              else
                discharge_grid[t,iState,jState] = discharge[t,iState,jState]+pv_prod[t]
              end
                  
            else                               #IDLING PHASE -> cal only sell power from PV to the grid (if any)
                
              charge[t,iState,jState] = 0
              discharge[t,iState,jState] = 0
              degradation[t,iState,jState] = 0

              charge_grid[t,iState,jState] = 0
              discharge_grid[t,iState,jState]=discharge[t,iState,jState]+pv_prod[t]

              if maximumGrid && discharge_grid[t,iState,jState]>= gridMax #se vi è un limite su potenza massima nella rete e la potenza da PV è maggiore della capacità di rete)
                discharge_grid[t,iState,jState] = gridMax
              else
                discharge_grid[t,iState,jState]=discharge[t,iState,jState]+pv_prod[t]
              end
      
            end

            val[t,iState,jState] = power_price[t]*NHoursStage*(discharge_grid[t,iState,jState]-charge_grid[t,iState,jState]) - degradation[t,iState,jState]*costBat*energy_Capacity - penaltySOC - penaltyExport + optimalValueStates[t+1,jState]      #/10E5
            gain[t,iState,jState] = power_price[t]*NHoursStage*(discharge_grid[t,iState,jState]-charge_grid[t,iState,jState])
            replacementCost[t,iState,jState] = degradation[t,iState,jState]*costBat*energy_Capacity      
            
            #= println(" State $iState")
            println(" State t:$iState, ",soc_start,"MWh - State t+1:$jState ",soc_final,"MWh - Charge[MW]: ",charge[t,iState,jState]," Discharge[MW]: ",discharge[t,iState,jState])
            println(" PV [MW]: ",pv_prod[t]," Charge from grid [MW]: ",charge_grid[t,iState,jState], " Discharge to grid[MW]: ",discharge_grid[t,iState,jState])
            println(" Value t+1 jState: ",optimalValueStates[t+1,jState])
            println(" Price [€/MWh]: ",power_price[t]," Value t [€]: ",val[t,iState,jState]," gain [€]: ",gain[t,iState,jState], " replacement cost [€]: ", replacementCost[t,iState,jState])
            println("") =#

          end # end jStates=1:5

        optimalValueStates[t,iState] = findmax(val[t,iState,:])[1]             # Trovo il massimo del Valore funzione obiettivo : transizioni + valore stato precedente 
        optimalfromState[t,iState] = findmax(val[t,iState,:])[2]               # Mi dice da quale stato al giorno precedente (o futuro) arrivo

        println("Optimal Val at stage t: $t and state x $iState: ",optimalValueStates[t,iState], " coming from state: ",optimalfromState[t,iState])
        println()

      end

    end   # end Stages

    # RACCOLGO I RISULTATI DEL PERCORSO MIGLIORE

    startingFrom = findmax(optimalValueStates[1,:])[2]                          # Mi inidica in quale stato per NStage=1 ho il massimo valore
    netOverallRevenues = 0 
    overallCost = 0

    for t=1:NStages   
      
      comingFrom = Int(optimalfromState[t,startingFrom])                        # Inidca da quale stato presedente sono arrivato
    
      optValue = findmax(optimalValueStates[t,startingFrom])[1]
      charge_bat = charge[t,startingFrom,comingFrom]
      dis_bat = discharge[t,startingFrom,comingFrom]
      char_grid = charge_grid[t,startingFrom,comingFrom]
      dis_grid = discharge_grid[t,startingFrom,comingFrom]
      pv = pv_prod[t]
      price = power_price[t]

      net_revenues = gain[t,startingFrom,comingFrom] - replacementCost[t,startingFrom,comingFrom]
      batCost = replacementCost[t,startingFrom,comingFrom]  

      overallCost = overallCost + batCost
      netOverallRevenues = netOverallRevenues + net_revenues

      push!(optimalPath,saveOptimalValues(t,optValue,price,startingFrom,comingFrom,seg,charge_bat,char_grid,dis_bat,dis_grid,pv,net_revenues,batCost))
      
      startingFrom=comingFrom

    end
    
    error = []
    for t=1:NStages 
      if optimalPath[t].charge> power_Capacity
        push!(error,optimalPath[t])
      end
    end

  end

  return Results_dp(
    charge,
    charge_grid,
    discharge,
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
    netOverallRevenues,
    error,
    modeFinal,
    costBat,
    gridMax,
   )

end


show(io::IO,x::optimalStage) = print(io,"Stage:",x.stage," -> opVal = ",x.optimalValue," price = ",x.price," curSOC = ",x.currentSOC,", nextSOC = ",x.nextSOC,", charge = ",x.charge, ", fromGrid = ",x.fromGrid,", discharge = ",x.discharge ,", toGrid = ",x.toGrid,", PV = ", x.PV ,", netRevenues = ",x.netRev, ", batCost = ",x.batteryCost)

function saveOptimalValues(stage::Int64,optimalValue::Float64,price::Float64, curSt::Int64, nextSt::Int64, seg::Any, charge::Float64,fromGrid::Float64,discharge::Float64, toGrid::Float64, PV::Float64, netRev::Float64 , batteryCost::Float64)
  optimalStage(stage, optimalValue, price, seg[curSt], seg[nextSt], charge, fromGrid, discharge,toGrid, PV,netRev, batteryCost)
end

#= ARROTONDANDO A DUE CIFRE SIGNIFICATIVE
function saveOptimalValues(stage::Int64,optimalValue::Float64,price::Float64, curSt::Int64, nextSt::Int64, seg::Any, charge::Float64,fromGrid::Float64,discharge::Float64, toGrid::Float64, PV::Float64, netRev::Float64 , batteryCost::Float64)
  optimalStage(stage, round(optimalValue;digits=2), round(price;digits=2), seg[curSt], seg[nextSt], round(charge;digits=2), round(fromGrid;digits=2), round(discharge;digits=2), round(toGrid;digits=2), PV, round(netRev;digits=2), round(batteryCost;digits=2))
end =#