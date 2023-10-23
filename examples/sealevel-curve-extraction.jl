# sealevel curve extraction

using Interpolations
using DataFrames
using CSV
using HTTP
using PyCall
dh = pyimport("datahugger")
url = "https://osf.io/5yka4/" #url of the folder
dh.get(url,"data/all-sealevel")

function sealevel()  
    #headers = Dict("Authorization" => "Bearer $token")
    #response = HTTP.get(url)# copy paste Sealelvel data from a database
    #data = split(String(response.body), "\n")
    #data_float64 = parse.(Float64, data)
    folder_path = "data/all-sealevel"
    files = readdir(folder_path)
    files = filter(x -> occursin(".txt", x), readdir(folder_path))

    for file in files
        file_path = joinpath(folder_path,file)
        name, _ = splitext(file)
        data = read(file_path, String)
        CSV.write("data/all-sealevel/$name.csv", DataFrame(sl = data))
    end

    #CSV.write("data/all-sealevel/$name.csv", DataFrame(sl=data))
end

#local_path = "data/all-sealevel/sl.csv"
sealevel() 
#=t = 1 : 2500
sl = sin.(2π .* t ./ 0.2)
#file = open("data/all-sealevel/sl1.csv","w")
CSV.write("data/all-sealevel/sl1.csv", DataFrame(sl=sl))
close(file)=#
