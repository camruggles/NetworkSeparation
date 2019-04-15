

function ReckfulRoaming(A,v,u,N,N2,RC)
    # use huda's code to compute the clustering coefficient
    for v in nodes
        CC_v = clustering coefficient for V
        RC_v = Int[]
        for u in N # get all of the neighbors of v
            # remove u from the neighbor list
            CC_rcV = new V clustering coefficient without u
            if CC_rcV > CC_v
                push!(RC_v, u)
            end
            reset neighbor nodes
        end
    for v in nodes
        for u in v neighbors
            for w in intersection of u and v # where does this node u come from?
                CC_uv = sum of clustering coefficients in the adjacency of u and v
                remove w from v's neighbors
                ? do we also remove
                CC_rvUV = sum of clustering coefficients in the adjacency of u and v ?without w?
                if CC_rvUV > CC_uv && Nv without w and Nu intersection is null && permit(v)
                    re add w to v's neighbors
                end
            end
            N_rrU = neighbors of v
        end
    end
end
