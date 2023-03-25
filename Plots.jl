# PLOTS
function plotPath(InputParameters::InputParam,ResultsDP::Results_dp,state_variables::states, priceYear)

    @unpack (NStages, NStates, NHoursStage, Big) = InputParameters
    @unpack (optimalPath,modeFinal) = ResultsDP;
    @unpack (seg) = state_variables

    optimalSOC=zeros(NStages+1)
    charge_battery = zeros(NStages)
    discharge_battery = zeros(NStages)
    charge_grid = zeros(NStages)
    discharge_grid = zeros(NStages)
    pv = zeros(NStages)

    for t=1:NStages
        optimalSOC[t]= optimalPath[t].currentSOC
        charge_battery[t] = optimalPath[t].charge
        discharge_battery[t] = optimalPath[t].discharge
        charge_grid[t] = optimalPath[t].fromGrid
        discharge_grid[t] = optimalPath[t].toGrid
        pv[t] = optimalPath[t].PV
    end
    optimalSOC[end]= optimalPath[end].nextSOC

    # PL0TTING OPTIMAL SOC FOR ALL PLANING PERIOD
    x=range(1, NStages+1, length=NStages+1)
    t1=plot(x,optimalSOC,label = "Optimal SOC", xlabel= "Stages",ylabel="SOC MWh", lc=:blue,size =(1500,600),yflip = false, show=true)
    display(t1)

    # PLOTTING CHARGE/DISCHARGE for BATTERY
    x1=range(1, NStages, length=NStages)
    t2=plot(x1,charge_battery-discharge_battery,label = "Charge-discharge", xlabel= "Stages",ylabel="Charge-discharge MW", lc=:blue,size =(1500,600),yflip = false, show=true)
    display(t2)
  
    #PLOTTING GRDI EXCHANGE
    x2=range(1, NStages, length=NStages)
    t3=plot(x2,charge_grid-discharge_grid,label = "Charge-discharge GRID", xlabel= "Stages",ylabel="Charge-discharge MW", lc=:blue,size =(1500,600),yflip = false, show=true)
    display(t3)

    savefig(t1, "Optimal SOC $modeFinal, $NStages stages $NStates states $priceYear.png")
    savefig(t2, "Charge-discharge $modeFinal, $NStages stages $NStates states $priceYear.png")
    savefig(t3, "Charge-discharge GRID $modeFinal, $NStages stages $NStates states $priceYear.png")

end
