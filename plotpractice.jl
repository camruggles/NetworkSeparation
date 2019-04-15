# plot number 8
# graph plots
# Now we can write our own graph plotting function!
@time begin
using LinearAlgebra
using SparseArrays
using Plots
using MatrixNetworks

# This function uses tricks we have learned so far to plot a graph
function graphplot(A,xy)
    f = plot(leg=false, axis = false)
    ei,ej,w = findnz(sparse(triu(A)))
    lx = [xy[ei,1]';xy[ej,1]';NaN*ones(1,length(ei))]
    ly = [xy[ei,2]';xy[ej,2]';NaN*ones(1,length(ei))]
    for i = 1:length(w)
        plot!(f,lx[:,i],ly[:,i],color = :black, linewidth = 1)
    end
    scatter!(f,xy[:,1],xy[:,2],color = :black)
    return f
end


using PyCall
pyimport_conda("igraph","python-igraph","conda-forge") # this will install if necessary
igraph = pyimport("igraph")
using MatrixNetworks
#A = load_matrix_network("cores_example")

function igraph_layout(A, layoutname::AbstractString="lgl")
    (ei,ej) = findnz(A)
    edgelist = [(ei[i]-1,ej[i]-1) for i = 1:length(ei)]
    nverts = size(A)
    G = igraph.Graph(nverts, edges=edgelist, directed=true)
    layoutname = "fr"
    xy = G[:layout](layoutname)
    xy = [ Float64(xy[i][j]) for i in 1:length(xy),  j in 1:length(xy[1])]
end

# include("verification.jl")
# Gt, Gert, ERt, q, q2, a2 = processFromFile()
# Gtp = Gt[q,q]
# Gertp = Gert[q,q]
# A = Gertp
n = 25
k = 3
k0 = 5
A = pa_graph(n, k, k0)
A = sparse_transpose(A)
println("1")
xy = igraph_layout(A)
println("2")
f = graphplot(A,xy)
println("3")
display(f)
end
println("4")
readline(stdin)

# draw out a graph and then compute the probabilities of what occurs
# implement the concept of raising the clustering coefficient

function preferential_attachment_edges!(
            nnew::Int,k::Int,edges::Vector{Tuple{Int,Int}},n0::Int)
    for iter=1:nnew
        println("\nedges1")
        @show edges
        i = n0+iter
        a = [rand(edges)[1] for j=1:k]
        println("\na")
        @show a
        newedges = unique(a)
        for v in newedges
            push!(edges, (i, v))
            push!(edges, (v, i))
        end
        println("\nedges2")
        @show edges
        @show length(edges)
        readline(stdin)
    end
    return edges
end

edges = []
for i = 1:5
    for v = 1:5
        if i != v
            push!(edges,(i,v))
        end
    end
end
edges = Array{Tuple{Int64, Int64}, 1}(edges)

preferential_attachment_edges!(20, 3, edges, 5)
