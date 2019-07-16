module Ca

# methods to export and make available
export foo, bar

foo(x::T, y::T) where T <: Real = x + y - 5
bar(z::Float64) = foo(sqrt(z), z)

# the initial correspondence analysis function
function ca(M::Array{Float64, 2}, k::Integer)
    O = M ./ sum(M)
    #row totals (row masses) by each of the column totals.
    E = sum(O, dims=2) .* sum(O, dims=1)
    Z = (O .- E) ./ sqrt.(E)
    U, S, V = svd(Z)
    U[:, k]
end

end # module
