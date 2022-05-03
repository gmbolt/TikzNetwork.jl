module TikzNetwork

using DataFrames, CSV, NetworkLayout, StatsBase

export get_tikz_input_files_muligraph, get_tikz_input_files_weighted

function relabel(edgelist::Vector{Tuple{Int,Int}})
    edgelist_relab = Tuple{Int,Int}[]
    vmap = Dict{Int,Int}()
    i = 1
    for edge in edgelist
        for v in edge
            if v ∉ keys(vmap)
                vmap[v] = i 
                i+=1
            end 
        end 
        push!(edgelist_relab, map(x->vmap[x], edge))
    end 
    rev_vmap = Dict(v=>k for (k,v) in vmap)
    return edgelist_relab, rev_vmap
end 

function relabel(edgelist::Vector{Tuple{String,String}})
    edgelist_relab = Tuple{Int,Int}[]
    vmap = Dict{String,Int}()
    i = 1
    for edge in edgelist
        for v in edge
            if v ∉ keys(vmap)
                vmap[v] = i 
                i+=1
            end 
        end 
        push!(edgelist_relab, map(x->vmap[x], edge))
    end 
    rev_vmap = Dict(v=>k for (k,v) in vmap)
    return edgelist_relab, rev_vmap
end 

function get_position(
    i::Int, 
    xlocs, ylocs,
    adjmat
    )
    
    ind = (adjmat[:,i] .> 0) .& ([j!=i for j in eachindex(xlocs)])
    if all(ylocs[ind] .< ylocs[i])
        return "above"
    elseif all(ylocs[ind] .> ylocs[i])
        return "below"
    else
        x_diff = mean(xlocs.-xlocs[i])
        y_diff = mean(ylocs.-ylocs[i])
        return angle_to_position(
            angle(complex(x_diff,y_diff) - complex(0,0))
        )
    end 
end

function get_label_angle(
    i::Int, 
    xlocs, ylocs,
    adjmat
    )
    
    ind = (adjmat[:,i] .> 0) .& ([j!=i for j in eachindex(xlocs)])
    if all(ylocs[ind] .< ylocs[i])
        return π/2
    elseif all(ylocs[ind] .> ylocs[i])
        return -π/2
    else
        x_diff = mean(xlocs.-xlocs[i])
        y_diff = mean(ylocs.-ylocs[i])
        return angle_to_position(
            angle(complex(0,0)-complex(x_diff,y_diff))
        )
    end 
end

function get_loop_angle(
    i::Int, 
    xlocs, ylocs,
    adjmat
    )
    
    ind = (adjmat[:,i] .> 0) .& ([j!=i for j in eachindex(xlocs)])
    if all(ylocs[ind] .< ylocs[i])
        return -90.0
    elseif all(ylocs[ind] .> ylocs[i])
        return 90.0
    else
        x_diff = mean(xlocs.-xlocs[i])
        y_diff = mean(ylocs.-ylocs[i])
        r = angle(complex(0,0)-complex(x_diff,y_diff))
        return rad2deg(r)
    end 
end

function angle_to_position(
    θ::Float64
    )

    if 0 < θ < π/8 
        return "left"
    elseif π/8 < θ < 3π/8
        return "below left"
    elseif 3π/8 < θ < 5π/8
        return "below"
    elseif 5π/8 < θ < 7π/8 
        return "below right"
    elseif 7π/8 < θ < π
        return "right"
    elseif -7π/8 < θ < -π
        return "right"
    elseif -5π/8 < θ < -7π/8
        return "above right"
    elseif -3π/8 < θ < -5π/8
        return "above"
    elseif -π/8 < θ < -3π/8
        return "above left"
    elseif 0 < θ < -π/8
        return "left"
    else 
        return "left"
    end 
end 

position_flip = Dict(
    "above" => "below",
    "above right" => "below left",
    "right" => "left",
    "below right" => "above left",
    "below" => "above",
    "below left" => "above right",
    "left" => "right",
    "above left" => "below right",
)

function get_tikz_input_files_muligraph(
    edges::Vector{Tuple{T,T}},
    vertices::Vector{T};
    filename_edges= "edges.csv",
    filename_vertices="vertices.csv",
    curves=true,
    bendstep=10,
    loopstep=0.5,
    vertexsize=0.35,
    directed=true,
    scale=3.0,
    number_edges=false,
    edge_color=:gray,
    edge_alpha = 1.0,
    edge_style="",
    vertex_color=:JuliaBlue,
    vertex_alpha=1.0,
    position="auto",
    loopshape=50,
    label_vertices=true
    ) where {T<:Union{Int,String}}

    if length(vertices)==length(edges)==0
        df_nodes = DataFrame()
        df_edges = DataFrame()
        CSV.write(filename_vertices, df_nodes)
        CSV.write(filename_edges, df_edges)
        return df_edges, df_nodes
    end 

    V = length(vertices)
    adjmat = zeros(Int, V, V)
    for (u,v) in edges
        i = findfirst(x->x==u, vertices)
        j = findfirst(x->x==v, vertices)
        adjmat[i,j] += 1 
    end 
    xlocs, ylocs = (Float64[],Float64[])
    locs = NetworkLayout.shell(adjmat)
    xlocs = [x[1] for x in locs] * scale
    ylocs = [x[2] for x in locs] * scale

    if V==1
        pos = ["above"]
    elseif position=="auto"
        pos = [get_position(i, xlocs, ylocs, adjmat) for i in eachindex(locs)]
    else 
        pos = [position for v in vertices]
    end 

    if typeof(vertex_color)==Symbol
        vcol = [string(vertex_color) for v in vertices]
    else 
        vcol = string.(vertex_color)
    end

    if length(vertex_alpha)==1
        valpha = [vertex_alpha[1] for v in vertices]
    else 
        valpha = vertex_alpha
    end

    df_nodes = DataFrame(
        id = vertices, 
        x=xlocs, 
        y=ylocs, 
        IdAsLabel=label_vertices,
        position=pos,
        size=[vertexsize for i in eachindex(vertices)],
        color=vcol,
        opacity=valpha
    )

    if length(edges)==0
        CSV.write(filename_vertices, df_nodes)
        df_edges = DataFrame()
        CSV.write(filename_edges, df_edges)
        return df_edges, df_nodes
    end 

    edges_seen = Dict{Tuple{T,T},Int}()
    if curves
        edge_curve_deg = []
        for e in edges 
            if e ∈ keys(edges_seen)
                edges_seen[e] += 1 
            else 
                edges_seen[e] = 1
            end 
            push!(edge_curve_deg, edges_seen[e]*bendstep)
        end 
    else 
        edge_curve_deg = [0 for e in edges]
    end 

    loopsize = Float64[]
    empty!(edges_seen)
    for e in edges
        if e[1]==e[2] #If self edge
            if e ∈ keys(edges_seen)
                edges_seen[e] += 1 
            else 
                edges_seen[e] = 1
            end 
            push!(loopsize, edges_seen[e]*loopstep)
        else
            push!(loopsize, loopstep)
        end 
    end 
    loop_pos = [get_loop_angle(findfirst(x->x==v, vertices), xlocs, ylocs, adjmat) for (u,v) in edges]

    if directed 
        is_directed = [true for e in edges]
    else 
        is_directed = [false for e in edges]
    end 

    if typeof(edge_color)==Symbol
        ecol = [string(edge_color) for e in edges]
    else 
        ecol = string.(edge_color)
    end

    if length(edge_alpha)==1
        ealpha = [edge_alpha[1] for e in edges]
    else 
        ealpha = edge_alpha
    end 

    if typeof(edge_style)==Symbol
        esty = [string(edge_style) for e in edges]
    else 
        esty = string.(edge_style)
    end 
    src = [e[1] for e in edges]
    dst = [e[2] for e in edges]

    df_edges = DataFrame(
        u=src, v=dst, 
        bend=edge_curve_deg, 
        Direct=is_directed,
        color=ecol,
        loopsize=loopsize,
        loopshape=loopshape,
        loopposition=loop_pos,
        opacity=edge_alpha,
        style=esty
    )

    if number_edges
        df_edges[:,"label"] = collect(eachindex(edges))
    end 

    CSV.write(filename_vertices, df_nodes)
    CSV.write(filename_edges, df_edges)

    return df_nodes, df_edges
end


function get_tikz_input_files_muligraph(
    edges::Vector{Tuple{T,T}};
    args...
    ) where {T<:Union{Int,String}}

    vertices = unique(vcat(collect.(edges)...))

    return get_tikz_input_files_muligraph(
        edges, vertices;
        args...
    )

end

function get_tikz_input_files_muligraph(
    vertices::Vector{T};
    filename_edges= "edges.csv",
    filename_vertices="vertices.csv",
    curves=true,
    bendstep=10,
    loopstep=0.5,
    vertexsize=0.35,
    directed=true,
    scale=3.0,
    number_edges=false,
    edge_color=:gray,
    position="auto"
    ) where {T<:Union{Int,String}}

    V = length(vertices)
    xlocs, ylocs = (Float64[],Float64[])
    locs = NetworkLayout.shell(zeros(Int, V, V))
    xlocs = [x[1] for x in locs] * scale
    ylocs = [x[2] for x in locs] * scale

    if V==1
        pos = ["above"]
    elseif position=="auto"
        pos = [get_position(i, xlocs, ylocs, adjmat) for i in eachindex(locs)]
    else 
        pos = [position for v in vertices]
    end 

    df_nodes = DataFrame(
        id = vertices, 
        x=xlocs, 
        y=ylocs, 
        IdAsLabel=true,
        position=pos,
        size=[vertexsize for i in eachindex(vertices)]
    )

    
    CSV.write(filename_vertices, df_nodes)
    df_edges = DataFrame()
    CSV.write(filename_edges, df_edges)
    return df_edges, df_nodes
    

end 


function get_tikz_input_files_weighted(
    edges::Vector{Tuple{T,T}},
    vertices::Vector{T};
    filename_edges= "edges.csv",
    filename_vertices="vertices.csv",
    bend=0,
    directed=true,
    vertexsize=0.35,
    scale=3.0,
    lwscale=2.0,
    edge_color=:gray,
    position="auto",
    labelvertices=true,
    loopsize=0.5, 
    loopshape=50
    ) where {T<:Union{Int,String}}

    V = length(vertices)
    adjmat = zeros(Int, V, V)
    for (u,v) in edges
        i = findfirst(x->x==u, vertices)
        j = findfirst(x->x==v, vertices)
        adjmat[i,j] += 1 
    end 

    xlocs, ylocs = (Float64[],Float64[])
    locs = NetworkLayout.shell(adjmat)
    xlocs = [x[1] for x in locs] * scale
    ylocs = [x[2] for x in locs] * scale

    if position=="auto"
        pos = [get_position(i, xlocs, ylocs, adjmat) for i in eachindex(locs)]
    else 
        pos = [position for v in vertices]
    end 

    edges_weighted = proportionmap(edges)

    if typeof(edge_color)==Symbol
        ecol = [string(edge_color) for e in keys(edges_weighted)]
    else 
        ecol = string.(edge_color)
    end

    if length(vertexsize)==1
        vsize = vertexsize[1]
    else 
        vsize = vertexsize
    end

    df_nodes = DataFrame(
        id = vertices, 
        x=xlocs, 
        y=ylocs, 
        IdAsLabel=labelvertices,
        position=pos,
        size=vsize
    )

    src = [e[1] for e in keys(edges_weighted)]
    dst = [e[2] for e in keys(edges_weighted)]

    weights = [edges_weighted[e] for e in keys(edges_weighted)]

    loop_pos = [get_loop_angle(findfirst(x->x==v, vertices), xlocs, ylocs, adjmat) for (u,v) in keys(edges_weighted)]

    df_edges = DataFrame(
        u=src, v=dst, 
        bend=bend, 
        Direct=[directed for e in keys(edges_weighted)],
        color=ecol,
        lw = weights * lwscale,
        opacity=0.8,
        loopposition=loop_pos,
        loopshape=loopshape,
        loopsize=loopsize
    )

    CSV.write(filename_vertices, df_nodes)
    CSV.write(filename_edges, df_edges)

    return df_nodes, df_edges
end





function get_tikz_input_files_weighted(
    edges::Vector{Tuple{T,T}};
    args...
    ) where {T<:Union{Int,String}}

    vertices = unique(vcat(collect.(edges)...))

    return get_tikz_input_files_weighted(
        edges, vertices;
        args...
    )

end



function get_tikz_input_files_weighted(
    edges::Vector{Vector{Tuple{T,T}}},
    vertices::Vector{T};
    filename_edges= "edges",
    filename_vertices="vertices",
    bend=0,
    directed=true,
    vertexsize=0.35,
    scale=3.0,
    lwscale=2.0,
    edge_color=:gray,
    position="auto",
    labelvertices=true,
    loopsize=0.5, 
    loopshape=50
    ) where {T<:Union{Int,String}}

    V = length(vertices)
    adjmat = zeros(Int, V, V)
    for (u,v) in vcat(edges...)
        i = findfirst(x->x==u, vertices)
        j = findfirst(x->x==v, vertices)
        adjmat[i,j] += 1 
    end 

    xlocs, ylocs = (Float64[],Float64[])
    locs = NetworkLayout.shell(adjmat)
    xlocs = [x[1] for x in locs] * scale
    ylocs = [x[2] for x in locs] * scale

    if position=="auto"
        pos = [get_position(i, xlocs, ylocs, adjmat) for i in eachindex(locs)]
    else 
        pos = [position for v in vertices]
    end 

    edges_weighted = [countmap(x) for x in edges]

    if typeof(edge_color)==Symbol
        ecol = [[string(edge_color) for e in keys(x)] for x in edges_weighted]
    else 
        ecol = [string.(x) for x in edge_color]
    end

    if length(vertexsize)==1
        vsize = vertexsize[1]
    else 
        vsize = vertexsize
    end

    df_nodes = DataFrame(
        id = vertices, 
        x=xlocs, 
        y=ylocs, 
        IdAsLabel=labelvertices,
        position=pos,
        size=vsize
    )

    src = [[e[1] for e in keys(x)] for x in edges_weighted]
    dst = [[e[2] for e in keys(x)] for x in edges_weighted]

    weights = [[x[e] for e in keys(x)] for x in edges_weighted]

    loop_pos = [
        [get_loop_angle(findfirst(x->x==v, vertices), xlocs, ylocs, adjmat) for (u,v) in keys(x)] for x in edges_weighted
    ]

    df_edges = [DataFrame(
        u=s, v=d, 
        bend=bend, 
        Direct=[directed for e in keys(x)],
        color=e,
        lw = w * lwscale,
        opacity=0.8,
        loopposition=l,
        loopshape=loopshape,
        loopsize=loopsize
    ) for (x,w,s,d,l,e) in zip(edges_weighted, weights, src, dst, loop_pos, ecol)]

    CSV.write(filename_vertices * ".csv", df_nodes)
    for (i,df) in enumerate(df_edges)
        CSV.write(filename_edges * "$i.csv", df)
    end 
    return df_nodes, df_edges
end

function get_tikz_input_files_weighted(
    edges::Vector{Vector{Tuple{T,T}}};
    args...
    ) where {T<:Union{Int,String}}

    vertices = unique(vcat(collect.(vcat(edges...))...))

    return get_tikz_input_files_weighted(
        edges, vertices;
        args...
    )

end

end