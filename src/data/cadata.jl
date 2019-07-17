export cadata

function cadata(name)
    basepath = joinpath(dirname(pathof(Ca)), "data")
    if name == "smoke"
        # get the data and convert to a matrix
        smoke = CSV.read(basepath*"/smoke.csv")
        names = convert(Array{String}, smoke[:, 1])
        dat = convert(Array{Int64, 2}, smoke[:, 2:end])
        return dat, names
    end

    if name == "author"
        #get the author data and convert to matrix
        author = CSV.read(basepath*"/author.csv")
        names = convert(Array{String}, author[:, 1])
        dat = convert(Array{Int64, 2}, author[:, 2:end])
        return dat, names
    end

end
