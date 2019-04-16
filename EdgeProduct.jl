using Plots

function computeEdgeProductROC(mixed, er, pa)
    n = size(mixed, 1)
    e = ones(n)

    nodeDegVector = mixed*e

    markerMatrix = mixed[:,:]
    for i = 1:n  # develop implementation that exploits the underlying sparse structure to go faster
        for j = 1:i
            if markerMatrix[i,j]==1 && er[i,j]==1
                markerMatrix[i,j] = 2
                markerMatrix[j,i] = 2
            end
        end
    end

    edgeVector = Array{Tuple{Float64,Tuple{Int64,Int64},Int64},1}()

    for i = 1:n
        for j = 1:i
            c = markerMatrix[i,j]
            if c != 0
                a = (i,j)
                b = nodeDegVector[i] * nodeDegVector[j]
                d = (b, a, c)
                push!(edgeVector, d)
            end
        end
    end

    sort!(edgeVector)

    tp::Int = 0
    tn::Int = 0
    fp::Int = 0
    fn::Int = 0

    sanityCheck = 0
    # t, tp, tn, fp, fn
    infoVector = Array{Tuple{Int, Int, Int, Int, Int},1}()

    for i in edgeVector
        a,b,c = i
        if c == 2
            tn += 1
        elseif c == 1
            fn += 1
        elseif c == 0
            println("there is an edge type of number id zero")
            println("this shouldn\'t happen, quitting program")
            exit(-1)
        end
        sanityCheck += 1
    end
    if sanityCheck != (tn + fn)
        println("sanity check failed")
        exit(-1)
    end
    push!(infoVector, (0, tp, tn, fp, fn))
    for t = 1:n
        a,b,c = edgeVector[t]
        if c==1
            tp += 1
            fn -= 1
        elseif c==2
            fp += 1
            tn -= 1
        end
        push!(infoVector, (t, tp, tn, fp, fn))
    end
    plotVector = Array{Tuple{Float64, Float64},1}()
    for entry in infoVector
        t, tp, tn, fp, fn = entry
        fpr = fp / (fp+tn)
        tpr = tp / (tp+fn)
        push!(plotVector, (fpr, tpr))
    end
    plot!(plotVector)
    readline("stdin")

end
