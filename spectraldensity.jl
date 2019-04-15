@time using Plots
using LinearAlgebra
using SparseArrays
using Arpack

@time include("verification.jl")

function plotSpectralDensity()
    A = rand(100,100)
    e = eigs(A, nev=100)
    h = histogram(e, nbins=100)
    display(h)
    println("Waiting")
    readline(stdin)
end

function plotGraphDensities()

    n = 1000
    k = 10
    k0 = 100
    p = 0.001

    @time Gert, Gt, ERt, q, q2, a2 = generateInputs(n,k,k0,p)
    q = q2 = a2 = 0 # get rid of some memory

    numEdges = sum(Gert.nzval)
    println("Mixed Edges: $numEdges")

    numEdges = sum(Gt.nzval)
    println("PA Edges: $numEdges")

    numEdges = sum(ERt.nzval)
    println("ER Edges: $numEdges")

    @time e1, ev = eigs(Gert, nev=n)
    @time e2, ev = eigs(Gt, nev=n)
    @time e3, ev = eigs(ERt, nev=n)

    data_labels1 = hcat("mixed graph")
    data_labels2 = hcat("preferential attachment")
    data_labels3 = hcat("erdos renyi")
    println("Histogram times")
    @time histogram(e1,nbins=100,fillalpha=0.3,labels=data_labels1)
    @time histogram!(e2,nbins=100,linecolor = :red,fillalpha=0.3,labels=data_labels2) # ,bar_position = :ss)
    @time h = histogram!(e3,nbins=100,linecolor = :green,fillalpha=0.3,labels=data_labels3) # ,bar_position = :ss)
    @time display(h)
    println("Waiting")
    readline(stdin)

end


@time plotGraphDensities()
