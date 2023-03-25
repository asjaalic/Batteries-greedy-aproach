# EXCEL SAVINGS
using DataFrames
using XLSX

function data_saving(InputParameters::InputParam,ResultsDP::Results_dp,Power_prices)

    @unpack (NStates,NStages,NHoursStage, Big) = InputParameters
    @unpack (optimalPath,modeFinal) = ResultsDP;

    hour=string(now())
    a=replace(hour,':'=> '-')

    nameF= "$modeFinal, $NStages stages, $NStates states, $a"

    #optimalState=zeros(NStages+1)
    optimalSOC=zeros(NStages+1);
    optimalChargeBattery=zeros(NStages+1);
    optimalDischargeBattery=zeros(NStages+1);
    optimalChargeGrid =zeros(NStages+1);
    optimalDischargeGrid=zeros(NStages+1);
    PV=zeros(NStages+1);
    price=zeros(NStages+1);
    degradationCost=zeros(NStages+1);
    netRevenues=zeros(NStages+1);

    for t=1:NStages
        optimalSOC[t]=optimalPath[t].currentSOC
        optimalChargeBattery[t]=optimalPath[t].charge
        optimalDischargeBattery[t]=optimalPath[t].discharge
        optimalChargeGrid[t]=optimalPath[t].fromGrid
        optimalDischargeGrid[t]=optimalPath[t].toGrid
        PV[t]=optimalPath[t].PV
        price[t]=optimalPath[t].price
        degradationCost[t]=optimalPath[t].batteryCost
        netRevenues[t]=optimalPath[t].netRev
    end
    optimalSOC[end]=optimalPath[end].nextSOC
    optimalChargeBattery[end] = 0
    optimalDischargeBattery[end] = 0
    optimalChargeGrid[end] = 0
    optimalDischargeGrid[end] = 0
    PV[end] = 0
    price[end] = Power_prices[NStages+1]
    degradationCost[end] = 0
    netRevenues[end] = 0

    table=DataFrame()
    table[!,"Stages"]= 1:1:(NStages+1)
    table[!,"Energy price €/MWh"] = price[:]
    table[!,"SOC MWh"] = optimalSOC[:]
    table[!,"Charge Battery MW"] = optimalChargeBattery[:]
    table[!,"Discharge Battery MW"]= optimalDischargeBattery[:]
    table[!,"PV production MW"] = PV[:]
    table[!,"Charge from Grid MW"] = optimalChargeGrid[:]
    table[!,"Discharge to Grid MW"] = optimalDischargeGrid[:]
    table[!,"Degradation Cost €"] = degradationCost[:]
    table[!,"Net Revenues €"] = netRevenues[:]

    XLSX.writetable("$nameF.xlsx", overwrite=true,
        results = (collect(DataFrames.eachcol(table)),DataFrames.names(table))
        )

end






