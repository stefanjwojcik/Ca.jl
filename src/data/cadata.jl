export cadata

function cadata(name)
    basepath = joinpath(dirname(pathof(Ca)), "data")
    if name == "smoke"
        # get the data and convert to a matrix
        dat = CSV.read(basepath*"/smoke.csv")
    end
    if name == "author"
        # get the data and convert to a matrix
        dat = CSV.read(basepath*"/author.csv")
    end
    if name == "haireye"
        # get the data and convert to a matrix
        dat = CSV.read(basepath*"/haireye.csv")
    end
    # create row and column names
    rnames = convert(Array{String}, dat[:, 1])
    cnames = names(dat)
    cnames = String.(cnames)[2:end]
    dat = convert(Array{Int64, 2}, dat[:, 2:end])
        return dat, rnames, cnames
end
