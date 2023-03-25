using CSV
using DataFrames
using DelimitedFiles

function read_csv(file, path)
    try
      println("Reading from: ", file)
      data = Matrix{Float64}(CSV.read(joinpath(path, file), DataFrame, header = false)) # €/MWh
      return data
    catch e
      println("Can't read ", file)
      throw(error())
    end
end

path = "C:\\Users\\Utente\\Desktop\\Batteries method2"
Power_prices = read_csv("prices_2020_8760.csv",path)
numData=length(Power_prices)

NewPrices=zeros(numData*2);

a=0
for i=1:(numData-1)
    NewPrices[i+a]=Power_prices[i]
    NewPrices[i+a+1] = (Power_prices[i]+Power_prices[i+1])/2
    a=a+1
end

NewPrices[end-1]=Power_prices[end]
NewPrices[end]= Power_prices[end]*0.9

writedlm("Prezzi mezz'ora 2020.csv", NewPrices, ' ')

