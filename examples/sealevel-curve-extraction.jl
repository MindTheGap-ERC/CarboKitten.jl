# sealevel curve extraction

using Interpolations
using DataFrames
using CSV
using HTTP
#import OpenScienceFramework as OSF

function sealevel(url,name)
    #=osf = OSF.Client(; token = alexa)  # put your OSF token here
    proj = OSF.project(osf; title="CarbonKitten")
    dir = OSF.directory(proj, "Analysis")
    subdir = OSF.directory(dir, "CarboKitten")
    subsubdir = OSF.directory(subdir, "OSF Storage")
    file = OSF.file(subsubdir, name)
    data = read(file)
    url = OSF.url(file)
    =#
    headers = Dict("Authorization" => "Bearer $token")
    response = HTTP.get(url,headers)# copy paste Sealelvel data from a database
    data = split(String(response.body), "\n")
    #data_float64 = parse.(Float64, data)
    CSV.write("data/all-sealevel/$name.csv", DataFrame(sl=data))
end

url = "https://osf.io/download/zntw8/"
#local_path = "data/all-sealevel/sl.csv"
name = "Auto000_Allo000_Stoch100V2"
token = "the token";
#localfile = "data/"
sealevel(url,name) 
#=t = 1 : 2500
sl = sin.(2π .* t ./ 0.2)
#file = open("data/all-sealevel/sl1.csv","w")
CSV.write("data/all-sealevel/sl1.csv", DataFrame(sl=sl))
close(file)=#
