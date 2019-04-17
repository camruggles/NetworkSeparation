
using LinearAlgebra
using SparseArrays
using MatrixNetworks
using Combinatorics
using Base
using Arpack
using DelimitedFiles
using Plots

#font = Plots.font("Helvetica", 10)
#plotly(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)
#plotly()

# show plot is the difference between 4 sec and 30 sec
showPlot = false # flag to determine if plotting will occur during this run
blocking = false
corrPrint = false


function testSeparation(PAS, ERS, PAO, ERO, n, q)

    # false positive, edge was added that wasn't there
    # false negative, edge wasn't added but was there

    # E = ERS-ERO
    #
    # e1 = count(i->(i==-1.0), E.nzval)
    # e2 = count(i->(i==1.0), E.nzval)
    # println("ERgraph: FP: $e2   FN: $e1")

    # count the number of edges in the original and separated graphs

    edgeCount1 = sum(PAS.nzval)
    edgeCount2 = sum(PAO.nzval)
    println("Number of edges in original vs new graphs")
    println("PA Edges Sep: $edgeCount1, orig: $edgeCount2")
    edgeCount1 = sum(ERS.nzval)
    edgeCount2 = sum(ERO.nzval)
    println("ER Edges Sep: $edgeCount1, orig: $edgeCount2")
    println("")

    # print out how many edges went to the wrong network during separation
    E = PAS-PAO
    e1 = count(i->(i==-1.0), E.nzval) # edges was not in separated but was in original
    e2 = count(i->(i==1.0), E.nzval) # edge was in separated but was not in original
    println("Incorrect separation placements")
    print("Edges from ER that went to PA: $e2\nEdges from PA that went to ER: $e1\n")
    println("")
    # no need to repeat for ERS and ERO, will have the same number as above
    # but flipped

    # plot the new degree distribution against the originals
    e = ones(n)
    a4 = PAS[q,q] * e # get pa distribution
    a5 = ERS[q,q] * e # get er distribution
    if showPlot
        if blocking
            println("Waiting") # a blocking operation to examine the original info
            readline(stdin)
        end

        data_labels = hcat("Separated pa distribution", "separated er distribution")
        g = Plots.plot!([a4 a5], labels=data_labels) # add degree distros to others
        display(g)

        if blocking
            println("Waiting") # blocking to prevent program from immediately terminating
            readline(stdin)
        end
    end

end


function generateInputs(n,k,k0,p)
    G = pa_graph(n, k, k0)
    Gt = sparse_transpose(G) # get matrix for pa graph

    q = degreeSort(Gt)  # index to permute the matrix so higher degree nodes have
    q2 = getOriginalPerm(q) # indices 1,2,.... etc.
    e = ones(n) # for row sum computes later on

    if showPlot
        C = Array{Int64}(Gt[q,q]) # get permuted matrix
        a1 = C*e # get degree distribution for original pa graph
    end

    # get matrix for ER graph
    ER = erdos_renyi_undirected(n,p)
    ERt = sparse_transpose(ER)
    # set overlapping edges as part of the PA graph
    @time for i in 1:n
        for j in 1:n
            if ERt[i,j] == 1 && ERt[i,j] == Gt[i,j]
                ERt[i,j] = 0
            end
        end
    end

    # get the mixed network and the permuted network and the degree distro
    Gert = Gt + ERt
    Gertp = Gert[q,q]
    a2 = Gertp * e

    # plot all of the degree distributions for the er, pa and mixed networks
    if showPlot
        a3 = ERt[q,q] * e # get the permuted distribution from the ER graph
        data_labels=hcat("original PA degree distribution", "original ER dist", "PA + ER degree distribution")
        f = Plots.plot([a1 a3 a2], labels=data_labels, title="n: $n | k: $k | k0: $k0 | p: $p")
        display(f)

    end
    #return the mixed and original networks, the index data for permutation,
    # and the mixed degree distribution
    dropzeros!(ERt)
    Gert, Gt, ERt, q, q2, a2

end

function randomClusterHypothesis(n,k,k0,p)
    @time Gert, Gt, ERt, q, q2, a2 = generateInputs(n,k,k0,p)
    Gtp = Gt[q,q]
    Gertp = Gert[q,q]

    @time x, x2 = plotClusterCoeffs(Gtp, Gertp, ERt, n, k, k0, p)

    x, x2


end

function colNormalized(M)
    A = 1.0*M
    n  =size(A,1)
    for j=1:n
        s = sum(A[:, j])
        A[:,j] /= s
    end
    @show A
    A

end

function main()

    # ((200*200)+(5000-200)*50) / 5000^2
    n = 5000
    k = 50
    k0 = 200
    p = 0.01

    # get mixed and original networks
    Gert, Gt, ERt, q, q2, a2 = generateInputs(n,k,k0,p)
    numEdges = length(ERt.nzval)

    # get permuted matrices, initialize empty er network
    Gertp = colNormalized(Gert[q,q])
    Gtp = colNormalized(Gt[q,q])
    b,c = eigs(Gertp,nev=1)
    if showPlot
        @show b
        plot(real.(c))
    end
    # @show c

    b,c = eigs(Gtp,nev=1)
    @show b
    if showPlot
        g = plot!(real.(c))
        display(g)
        readline(stdin)
    end
    # @show c


end

function getClusterCoeff(A)

    clustercoeffs(A,false,true)
end

function generateToFile()
    n = 5000
    k = 50
    k0 = 200
    p = 0.01

    # get mixed and original networks
    Gert, Gt, ERt, q, q2, a2 = generateInputs(n,k,k0,p)
    writedlm("BA1.dlm", Gt)
    writedlm("ER1.dlm", ERt)
    writedlm("Mixed1.dlm", Gert)
    B = [q q2 a2]
    writedlm("graph1info.dlm", B)

end

function processFromFile()
    Gt = sparse(Array{Float64}(readdlm("BA1.dlm")))
    Gert = sparse(Array{Float64}(readdlm("Mixed1.dlm")))
    ERt = sparse(Array{Float64}(readdlm("ER1.dlm")))
    B = readdlm("graph1info.dlm")
    q = B[:, 1]
    q2 = B[:, 2]
    a2 = B[:, 3]

    q = Int.(q)
    q2 = Int.(q2)

    Gt, Gert, ERt, q, q2, a2

end


function clusterHypothesis2()
    Gt, Gert, ERt, q, q2, a2 = processFromFile()
    Gtp = Gt[q,q]
    Gertp = Gert[q,q]

    A = Gtp[205:5000, 205:5000]
    B = Gertp[205:5000, 205:5000]

    x, x2 = plotClusterCoeffs(A, B)
    readline(stdin)

    x, x2


end

function clusterHypothesis3()
    Gt, Gert, ERt, q, q2, a2 = processFromFile()
    Gtp = Gt[q,q]
    Gertp = Gert[q,q]

    C = Gtp[1:200, :] # 200x5000 = 1x5000
    e = ones(1, 200)
    Ccolsums = e*C
    for i in 1:5000
        if Ccolsums[i] > 0
            Ccolsums[i] = 1
        end
    end
    Ccolsums = Ccolsums[:, 201:5000]

    D = Gertp[1:200, :] # 200x5000 = 1x5000
    # e = ones(1,200)
    Dcolsums = e*D
    for i in 1:5000
        if Dcolsums[i] > 0
            Dcolsums[i] = 1
        end
    end
    Dcolsums = Dcolsums[:, 201:5000]


    A = Gtp[201:5000, 201:5000]
    B = Gertp[201:5000, 201:5000]
    A = sparse([0 Ccolsums; Ccolsums' A])
    B = sparse([0 Dcolsums; Dcolsums' B])

    x, x2 = plotClusterCoeffs(A, B)
    readline(stdin)

    x, x2
end

function clusterHypothesis()
    Gt, Gert, ERt, q, q2, a2 = processFromFile()
    Gtp = Gt[q,q]
    Gertp = Gert[q,q]

    x, x2 = plotClusterCoeffs(Gtp, Gertp)

    x, x2


end


function plotClusterCoeffs(Gt, Gert)
    # ? remove the largest node cluster
    x = getClusterCoeff(Gt)
    x2 = getClusterCoeff(Gert)
    n = length(x)
    # x = x[2:n]
    # x2 = x2[2:n]
    f = plot([x x2], labels=hcat("unmixed", "mixed"))
    display(f)
    readline(stdin)

    x, x2
end

function plotClusterCoeffs(Gt, Gert, ERt, n, k, k0, p)
    # ? remove the largest node cluster
    x = getClusterCoeff(Gt)
    numEdges = sum(Gt.nzval)/2
    numRandomEdges = sum(ERt.nzval)/2
    x2 = getClusterCoeff(Gert)
    n = length(x)
    # x = x[2:n]
    # x2 = x2[2:n]
    f = plot([x x2], labels=hcat("unmixed", "mixed"), title="n: $n | k: $k | k0: $k0 | p: $p \n edges: $numEdges | randomEdges: $numRandomEdges")
    display(f)

    #readline(stdin)

    x, x2
end

function oldMain()
        # causes overflow
        # n = 100
        # k = 5
        # k0 = 20
        # p = 0.01


            # n = 100
            # k = 2
            # k0 = 10
            # p = 0.01

    n = 500
    k = 8
    k0 = 30
    p = 0.01

    # get mixed and original networks
    Gert, Gt, ERt, q, q2, a2 = generateInputs(n,k,k0,p)
    numEdges = length(ERt.nzval)

    # get permuted matrices, initialize empty er network
    Gertp = Gert[q,q]
    newErt = zeros(Int, n,n)

    # Given a matrix, compute the full correlation matrix
    # and also find the edge of the smallest probability
    # and separate it
    lowerBoundCorrelation = 0.0
    e = ones(n)
    println("Time to run all iterations")
    for z = 1:(numEdges/2)
        # initialized information for finding minimum degree correlation
        minIndex = (0,0)
        minValue = 1
        # nodeCorrelationMatrix = zeros(n,n)

        # for all edges calculate the correlation
        for i = 1:n
            for j = 1:n

                if Gertp[i,j] == 0
                    continue
                end

                # get the node degree correlation
                l = a2[j]
                l2 = a2[i]
                m = k

                val = BANodeCorrelation(l, l2, m)
                # nodeCorrelationMatrix[i,j] = val
                #@show val

                # check if it has the smallest node degree correlation
                if val < minValue
                    minValue = val
                    minIndex = (i,j)
                end

            end
        end
        # print for informative purposes
        lowerBoundCorrelation = minValue
        if corrPrint
            @show lowerBoundCorrelation
        end

        #@show minIndex
        i,j = minIndex

        # defensive coding, make sure it didn't get a zero edge
        if Gertp[i,j] == 0
            println("Critical error, node correlation positive when there is no edge")
            return -1
        end

        # remove the edge from the mixed network and put it in the er network
        Gertp[i,j] = 0
        Gertp[j,i] = 0
        newErt[i,j] = 1
        newErt[j,i] = 1

        # get a new degree distribution for calculating node correlations
        a2 = Gertp * e
    end
    println("")


    println("Smallest value of node degree correlations")
    @show lowerBoundCorrelation
    println("")

    # get original matrix permutations
    Gertp = Gertp[q2, q2]
    newErt = sparse(newErt[q2, q2])

    # get the number of edges in the new pa matrix
    println("Number of nodes in new matrix")
    @show sum(Gertp*e)
    println("")

    # show norm differences
    println("Norms")
    @show norm(Gert-Gertp)
    @show norm(newErt-ERt)
    println("")

    # show other analysis
    testSeparation(Gertp, newErt, Gt, ERt, n, q)

# why are NOTE AND BUG highlighted?
# IMPORTANT NOTE OF BUG FOUND EARLIER
#final gertp has the same number of edges as the original gert
# this is not supposed to happen, as edges were supposed to be separated
# turns out when I switched a one to a zero, nzval was still storing the zero
# do sum(nzval) not length(nzval)


end

# function to compute node degree correlations in a ba model
function BANodeCorrelation(l2, k2, m)
    # https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model
    l = Int(l2)
    k = Int(k2)
    # P (l | k) with m connection upon entry
    # In other words, if we select a node with degree k, and then select one neighbor randomly
    # the probability that this randomly selected neighbor will have degree L is given
    # by the function.
    a::BigInt = Base.binomial(BigInt(2*m+2), BigInt(m+1))
    b::BigInt = Base.binomial(BigInt(k+l-2m), BigInt(l-m))
    c::BigInt = Base.binomial(BigInt(k+l+2), BigInt(l+1))
    ret = m*(k+2) / (k*l*(l+1)) * (1 - a*b/c )
    if ret > 1.0
        @show ret
        println("Probability greater than one, you broke math")
        println("Program terminated early")
        Base.exit()
    end
    ret
end

# get the inverse of a matrix permutation
function getOriginalPerm(p)
    n = length(p)
    a = zeros(Int32, n)
    for i = 1:n
        x = p[i]
        a[x] = Int(i)
    end
    a
end
#=
example permuation:
[ 3 5 4 1 2]

example of unpermutation
[ 4 5 1 3 2  ]

unperm[1] = i where perm[i] = 1

=#


# get a permutation of a matrix
# higher degrees have lower indices
# if i < j, deg(A[i]) > deg(A[j])
function degreeSort(Gt)
    # getting the distribution of degrees in an array
    n = size(Gt,1)
    e = ones(n)
    degDistro = Gt*e
    nodeLabels = ones(n)

    # create a list of tuples of the degree count with the original node label
    tupleList = []
    for i = 1:n
        nodeLabels[i] = i
        push!(tupleList, (degDistro[i], nodeLabels[i]))
    end

    # sort the tuple List so higher degrees come first
    sort!(tupleList, rev=true)

    # create a relabelling s.t. 1 gets the high degree distribution
    p = vec(zeros(Int64, n))
    for i = 1:n
        a,b = tupleList[i]
        p[i] = Int(b)
    end

    p


end

# function to directly create a sparse matrix after degree sorting
# do this without using A = A[q,q]
# not tested
function degreeSortedSparse(n, nnz, sortedNodes)
    colptr = zeros(n)
    rowvals = zeros(nnz)
    nzval = zeros(nnz)
    k=1
    ncolvals=1
    for i in 1:n
        colptr[i] = ncolvals
        for j in 1:n
            si = sortedNodes[i]
            sj = sortedNodes[j]
            if G[si,sj] == 0
                continue
            else
                # push!(rowvals, j)
                rowvals[k] = j
                nzval[k] = G[si,sj]
                k += 1
                # push!(nzval, j)
            end
        end
    end
end

# convert a matrix of booleans to ones and zeros
# there is likely a julia function that will do this automatically
# I just don't know it
function BoolToInt(A)
    B = zeros(Int64, n, n)
    for i in 1:n
        for j in 1:n
           if A[i,j]
               B[i,j] = 1
           end
       end
    end
    B
end

# create an integer matrix from a boolean matrix
function intMat(n)
    A = rand(Bool, n, n)
    B = zeros(Int64, n, n)
    for i in 1:n
        for j in 1:n
           if A[i,j]
               B[i,j] = 1
           end
       end
    end
    A=B
    A
end

# O(n!) sorting algorithm for S's and G's
function badSort(a)
    b = permutations(a)
    @show length(b)
    for arr in b
        n = length(arr)
        sorted = true
        for i = 1:n-1
            if arr[i] > arr[i+1]
                sorted = false
                break
            end
        end
        if sorted
            return copy(arr)
        end
    end
    zeros(length(a))
end

#===========================================#
#clusterHypothesis()
