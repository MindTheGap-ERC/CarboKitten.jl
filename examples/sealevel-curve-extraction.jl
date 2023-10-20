# sealevel curve extraction
#module Import_SL

using Interpolations
using DataFrames
using CSV
using HTTP
#import OpenScienceFramework as OSF

function sealevel(url)
    #=osf = OSF.Client(; token = alexa)  # put your OSF token here
    proj = OSF.project(osf; title="MIND THE GAP")
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
    CSV.write("data/all-sealevel/sl.csv", DataFrame(sl=data))
end

#url = "https://raw.githubusercontent.com/xyl96/sealevel_curve/main/stoch1" # extracting from an open url
url = "https://osf.io/download/b956h/"
#local_path = "data/all-sealevel/sl.csv"
#name = "Auto000_Allo000_Stoch100V2.txt"
token = "ohXRgQVgTj0bSR8hA9mItYdWVnK2m4gyRDjFFrfJDkHFlmt7PPVNFNFj0xtQUThLMSQUIJ";
#localfile = "data/"
sealevel(url) 
#=t = 1 : 2500
sl = sin.(2π .* t ./ 0.2)
#file = open("data/all-sealevel/sl1.csv","w")
CSV.write("data/all-sealevel/sl1.csv", DataFrame(sl=sl))
close(file)=#
