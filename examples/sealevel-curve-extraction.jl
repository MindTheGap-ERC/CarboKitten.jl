# sealevel curve extraction
#module Import_SL

using Interpolations
using DataFrames
using CSV
using HTTP

function sealevel(url)
    #headers = Dict("Authorization" => "Bearer $token")
    response = HTTP.get(url)# copy paste Sealelvel data from a database
    data = split(String(response.body), "\n")
    data_float64 = parse.(Float64, data)
    CSV.write("data/all-sealevel/sl.csv", DataFrame(sl=data_float64))

end

url = "https://raw.githubusercontent.com/xyl96/sealevel_curve/main/stoch1" # extracting from an open url
#local_path = "data/all-sealevel/sl.csv"
#token = "ohXRgQVgTj0bSR8hA9mItYdWVnK2m4gyRDjFFrfJDkHFlmt7PPVNFNFj0xtQUThLMSQUIJ"
#localfile = "data/"
sealevel(url) 
#=t = 1 : 2500
sl = sin.(2π .* t ./ 0.2)
#file = open("data/all-sealevel/sl1.csv","w")
CSV.write("data/all-sealevel/sl1.csv", DataFrame(sl=sl))
close(file)=#
