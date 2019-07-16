
using BenchmarkTools;
using PyCall;
using CSV;
using LinearAlgebra;
using RCall;
using DataFrames

df = CSV.read("/Users/stefanwojcik/Documents/Julia Stuff/bots_correspondence_large_example.csv") |> DataFrame

df = convert(Array{Float64, 2}, df[:, 2:end])

#convert(Matrix{Float64}, jtab);

mutable struct correspondence
    M::Array{Float64,2}
    mat_total::Float64
    O::Array{Float64,2}
    E::Array{Float64,2}
    Z::Array{Float64,2}
    U::Array{Float64,2}
    S::Array{Float64,1}
    V::Array{Float64,2}
end

function correspondence(sizes)
    M = randn(sizes)
    mat_total = 0
    O = randn(sizes)
    E = randn(sizes)
    Z = randn(sizes)
    U = randn(sizes)
    S = randn(sizes[1])
    V = randn(sizes)
    correspondence(M, mat_total, O, E, Z, U, S, V)
end

function ca_in_julia(ca::correspondence, input_mat::Array{Float64, 2}, k::Int)
    ca.M[:] = input_mat
    ca.mat_total = sum(ca.M)
    ca.O[:] = ca.M ./ ca.mat_total
    #row totals (row masses) by each of the column totals.
    ca.E[:] = sum(ca.O, dims=2) .* sum(ca.O, dims=1)
    ca.Z[:] = (ca.O .- ca.E) ./ sqrt.(ca.E)
    ca.U[:], ca.S[:], ca.V[:] = svd(ca.Z)
    ca.U[:, k]
end

ca_ex = correspondence((10,10))
ca_in_julia(ca_ex, randn(10,10))


function ca_in_juliaV1(M::Matrix{Float64}, k::Integer)
    O = M ./ sum(M)
    #row totals (row masses) by each of the column totals.
    E = sum(O, dims=2) .* sum(O, dims=1)
    Z = (O .- E) ./ sqrt.(E)
    U, S, V = svd(Z)
    out = U[:, k]
end


# THE R CODE write a correspondence analysis function: m is a person-document matrix, k is number of dimensions
R"""
library(ca)
new_table = read.csv('/users/stefanwojcik/Documents/Julia\ Stuff/bots_correspondence_large_example.csv')
tab = new_table[, 2:ncol(new_table)]
rownames(tab) = new_table$user_name
""";

jtab = @rget tab;
jtab = convert(Matrix{Float64}, jtab);
M = jtab
