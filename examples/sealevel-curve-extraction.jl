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

    folder_path = "data/all-sealevel"
    files = readdir(folder_path)
    files = filter(x -> occursin(".txt", x), readdir(folder_path))

    for file in files
        file_path = joinpath(folder_path,file)
        name, _ = splitext(file)
        data = split(read(file_path, String),"\n")
        CSV.write("data/all-sealevel/$name.csv", DataFrame(sl = data))
    end

end

sealevel() 

