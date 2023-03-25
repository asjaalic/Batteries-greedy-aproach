using Rainflow
using Plots

# FIND ALL POSSIBLE TRAJECTORIES - combinations k^n 

power = 80                                                                                  # MW
NStages = 9                                                                                 # considero NStages come numero di stadi ll'interno del problema               
NStates = 8                                                                                 # per ogni statio ho NStates possibili stati
HoursPerStages = 24/(NStages-1)                                                             # hours per stage
@time newVector=collect.(Iterators.product(ntuple(_ -> 1:NStates, NStages)...))[:]          # vettore che contiene tutte le possibile combinazioni Nstates^Nstages
#Vector2=collect.(Iterators.product(ntuple(_ -> 1:4, 13)...))[:]

SOC =zeros(NStages,NStates)
SOCperc = zeros(NStages,NStates)
prices= rand(NStages-1)*100                                                                 #  presso al €/MWh
cost = rand(1)*100                                                                          #  genera un numero random per il costo della batteria (assumiamo costante)

for t=1:NStages
    for iState=1:NStates
        a= round(power- power*(1/(NStates-1)*(iState-1)); digits=2)
        SOC[t,iState] = a
        b= round(100*a/power; digits=2)
        SOCperc[t,iState]=b
        println("Stage $t , state $iState , $a MW , in percentage $b %")
    end
end

combinations=length(newVector)
num=rand(1:combinations)
plot(num, title= "States of power", label= "Trajectory $num",color=:blue,legend=:topleft)
plot!(twinx(),prices, label = "Energy price €/MWh",color=:red)
extremes,t=sort_peaks_NEW(num)
scatter!(t,extremes)

powerExchanged = zeros(combinations,NStages-1)
states_power =zeros(combinations,NStages)
gain = zeros(combinations,NStages-1)
totGain=zeros(combinations)

for j=1:combinations
    for t=1:NStages-1
        state_t= newVector[j][t]
        state_next_t= newVector[j][t+1]
        powerExchanged[j,t]=SOC[t+1,state_next_t]-SOC[t,state_t]
        gain[j,t]=-powerExchanged[j,t]*prices[t]*HoursPerStages
    end
    for t=1:NStages
        state_t= newVector[j][t]
        states_power[j,t]=SOC[t,state_t]
    end
    totGain[j]=sum(gain[j,:])
end


# COMPUTE THE RAINFLOW ALGORITHM and FIND EQUIVALENT CYCLES FOR EVERY TRAJECTORY

totCycles=zeros(a)

signal=states_power[14300,:]
extremes,t=sort_peaks_NEW(signal)
plot(signal)
scatter!(t,extremes)
cycle=count_cycles(extremes,t)
tot_cycles=sum_cycles(cycle,1,1)


for i=1:a
    segnali=states_power[i,:]
    estremi, t = sort_peaks_NEW(segnali)
    cycles=count_cycles(estremi,t)
    lunghezza=length(cycles)
    for j=1:lunghezza
        totCycles[i]=totCycles[i]+cycles[j].count
    end
end

#  Find trajectory that minimizes the overall costs considering battery replacement
finalCosts=zeros(a)

for i=1:a
    finalCosts[i]=totGain[i].-totCycles[i]*cost[1]
end

# FIND TRAJECTORY WITH MAXIMUM GAIN
findmax(finalCosts)
s=findmax(finalCosts)[2]
signal=states_power[s,:]
extremes,t=sort_peaks_NEW(signal)
plot(signal,label="States of power",legend=:topleft)
scatter!(t,extremes)
cycle=count_cycles(extremes,t)
powerExchanged[s,:]
prices
plot!(twinx(),prices,color=:red,xticks=:none,label="price [€/MWh]")
gain[s,:]

# FIND TRAJECTORY WITH MINIMUM GAIN
findmin(finalCosts)
p=findmin(finalCosts)[2]
signal_min=states_power[p,:]
extremes,t=sort_peaks_NEW(signal_min)
plot(signal_min)
scatter!(t,extremes)
cycle=count_cycles(extremes,t)


#=
random_trajectory = rand(1:a)
signal=states_power[random_trajectory,:]
extremum, t = sort_peaks(signal)
plot(signal, xticks=1:NStages)
scatter!(t,extremum) # plots extremes
cycles = count_cycles(extremum, t)
lunghezza=length(cycles)
totC=0
    for j=1:lunghezza
        totC=totC+cycles[j].count
    end
totC
=#


function sort_peaks_NEW(signal::AbstractArray{Float64,1}, dt=collect(1.:length(signal)))
    slope = diff(signal)
    # Determines if the point is local extremum
    is_extremum = vcat(true, (slope[1:end-1].*slope[2:end]).<0., true)
    return signal[is_extremum] , dt[is_extremum]
end

signal = 10*randn(100); # Gennerates some random data
extremum, t = sort_peaks(signal) # Sorts the signal to local extrema's, could optional take a time vector
plot(signal)
scatter!(t,extremum) # plots extrema's
cycles = count_cycles(extremum, t) # find all the cycles in the data
bins = sum_cycles(cycles, 10, 1) # Sums the cycles together dependant on intervals, here there is 10 range interval and one mean interval
bar(bins)