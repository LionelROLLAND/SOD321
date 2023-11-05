using JuMP, Gurobi, Graphs, GraphRecipes, Plots, Colors, GraphPlot, FixedPointNumbers, MetaGraphs

ENV["GRB_LICENSE_FILE"] = "/opt/gurobi1003/linux64/gurobi.lic"
include(joinpath("src", "read.jl"))


#### Structures de données utiles


struct GraphInfo
    neighb::Vector{Set{Int64}}
    arc_values::Any
end

struct ConstraintData
    n::Int64
    d::Int64
    f::Int64
    Amin::Int64
    Nr::Int64
    regions_to_nodes::Dict{Int64,Vector{Int64}}
    D::Matrix{Int64}
end

struct PbData
    arcs::Set{Tuple{Int64,Int64}}
    neighb::Vector{Set{Int64}}
    rev_neighb::Vector{Set{Int64}}
    c_data::ConstraintData
end

mutable struct TarjanState
    fifo::Vector{Int64}
    is_in_fifo::Vector{Bool}
    first_seen_time::Vector{Int64}
    lowest_accessible::Int64
    time::Int64
    components::Vector{Vector{Int64}}
end

struct Path
    nodes::Vector{Int64}
    length::Int64
end


#### Fonctions utiles


function create_pb_data(instance::String)::PbData
    n, d, f, Amin, Nr, R, regions, coords, D = readInstance(instance)

    neighb = Vector{Set{Int64}}(undef, n)
    rev_neighb = Vector{Set{Int64}}(undef, n)
    arcs = Set{Tuple{Int64,Int64}}()

    for i in 1:n
        neighb[i] = Set{Int64}()
        rev_neighb[i] = Set{Int64}()
    end

    for i in 1:n
        for j in 1:n
            if D[i, j] <= R && j != i && i != f && j != d
                push!(neighb[i], j)
                push!(rev_neighb[j], i)
                push!(arcs, (i, j))
            end
        end
    end

    regions_to_delete = Vector{Int64}(undef, 0)

    for (r, nodes) in pairs(regions)
        if d in nodes || f in nodes
            push!(regions_to_delete, r)
        end
    end
    for r in regions_to_delete
        delete!(regions, r)
    end
    Nr = length(regions)

    c_data = ConstraintData(
        n,
        d,
        f,
        Amin,
        Nr,
        regions,
        D,
    )

    res = PbData(
        arcs,
        neighb,
        rev_neighb,
        c_data,
    )
    return res
end

function flow_creation(; node::Int64, c_data::ConstraintData)
    if node == c_data.d
        return 1
    elseif node == c_data.f
        return -1
    else
        return 0
    end
end

function makepath(chosen_arcs::Set{Tuple{Int64,Int64}})::Vector{Int64}
    arc_dict = Dict{Int64,Int64}()
    for (i, j) in chosen_arcs
        if i in keys(arc_dict)
            println("Error : The set of chosen arcs doesn't define an elementary path.")
            return [-1]
        end
        arc_dict[i] = j
    end
    s = only(setdiff(keys(arc_dict), values(arc_dict)))
    current_end = s
    res = Vector{Int64}(undef, 0)
    while current_end != -1
        push!(res, current_end)
        current_end = pop!(arc_dict, current_end, -1)
    end
    if !isempty(arc_dict)
        println("Error : The set of chosen arcs doesn't define an elementary path.")
        return [-1]
    end
    return res
end



#### Resolution avec nombre de contraintes polynomial


function poly_solve(pb_data::PbData)::Path
    neighb = pb_data.neighb
    rev_neighb = pb_data.rev_neighb
    arcs = pb_data.arcs
    c_data = pb_data.c_data
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_time_limit_sec(model, 120.0)
    @variable(model, a[arcs], Bin)
    @variable(model, t[setdiff(1:c_data.n, (c_data.d, c_data.f))])
    @constraint(
        model,
        flow_cons[i in 1:c_data.n],
        sum(a[(i, j)] for j in neighb[i])
        -
        sum(a[(j, i)] for j in rev_neighb[i])
        ==
        flow_creation(node=i, c_data=c_data),
    )  # conservation du flot
    @constraint(model, bornitude[i in 1:c_data.n], sum(a[(i, j)] for j in neighb[i]) <= 1)  # Bornitude
    @constraint(
        model,
        exit_region[r in keys(c_data.regions_to_nodes)],
        sum(
            a[(i, j)] for i in c_data.regions_to_nodes[r]
            for j in setdiff(neighb[i], c_data.regions_to_nodes[r])
        ) >= 1,
    )  # Sortie de chacune des regions
    @constraint(model, min_visits, sum(a) >= c_data.Amin - 1)  # Visite au moins Amin aerodromes
    @constraint(
        model,
        elementarite[(i, j) in filter(((i, j),) -> isdisjoint((c_data.d, c_data.f), (i, j)), arcs)],
        t[j] >= t[i] + 1 + (c_data.n - 1) * (a[(i, j)] - 1),
    )  # Contrainte pour assurer l'elementarite (MTZ)
    @objective(model, Min, sum(c_data.D[i, j] * a[(i, j)] for (i, j) in arcs))
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    @assert primal_status(model) == FEASIBLE_POINT
    chosen_arcs = Set{Tuple{Int64,Int64}}((i, j) for (i, j) in arcs if value(a[(i, j)]) > 0.5)
    path_length = sum(c_data.D[i, j] for (i, j) in chosen_arcs)
    path = Path(makepath(chosen_arcs), path_length)
end



#### Sous-problème pour les contraintes GCS



epsilon = 0.000_001

function separation_pb(; pb_data::PbData, tilde_a::Dict{Tuple{Int64,Int64},Float64})::Tuple{Float64,Set{Int64}}
    neighb = pb_data.neighb
    rev_neighb = pb_data.rev_neighb
    arcs = pb_data.arcs
    c_data = pb_data.c_data
    chosen_arcs = Set{Tuple{Int64,Int64}}(
        (i, j) for (i, j) in arcs
        if value(tilde_a[(i, j)]) > epsilon
        &&
        isdisjoint((c_data.d, c_data.f), (i, j))
    )
    tilde_big_n = Set{Int64}(l for (i, j) in chosen_arcs for l in (i, j))
    sep_model = Model(Gurobi.Optimizer)
    set_silent(sep_model)
    @variable(sep_model, b[i in tilde_big_n], Bin)
    @variable(sep_model, u[(i, j) in chosen_arcs], Bin)
    @variable(sep_model, k[i in tilde_big_n], Bin)
    @constraint(sep_model, [i in tilde_big_n], k[i] <= b[i])  # Def de k
    @constraint(sep_model, [(i, j) in chosen_arcs], u[(i, j)] <= b[i])  # Def de u (1)
    @constraint(sep_model, [(i, j) in chosen_arcs], u[(i, j)] <= b[j])  # Def de u (2)
    @constraint(sep_model, sum(k) == 1)  # Un seul k choisi
    @constraint(sep_model, sum(b) >= 2)  # |S| >= 2
    @objective(
        sep_model,
        Max,
        sum(tilde_a[(i, j)] * u[(i, j)] for (i, j) in chosen_arcs)
        -
        sum(sum(tilde_a[(i, j)] for j in neighb[i]) * (b[i] - k[i]) for i in tilde_big_n)
    )
    optimize!(sep_model)
    @assert termination_status(sep_model) == OPTIMAL
    @assert primal_status(sep_model) == FEASIBLE_POINT
    return objective_value(sep_model), Set{Int64}(i for i in tilde_big_n if value(b[i]) > 0.5)
end



#### Resolution total du problème avec contraintes GCS + séparation PLNE



function expo_solve(pb_data::PbData)::Path
    neighb = pb_data.neighb
    rev_neighb = pb_data.rev_neighb
    arcs = pb_data.arcs
    c_data = pb_data.c_data
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    @variable(model, 0 <= a[arcs] <= 1, Bin)
    @constraint(
        model,
        flow_cons[i in 1:c_data.n],
        sum(a[(i, j)] for j in neighb[i])
        -
        sum(a[(j, i)] for j in rev_neighb[i])
        ==
        flow_creation(node=i, c_data=c_data),
    )  # conservation du flot
    @constraint(model, bornitude[i in 1:c_data.n], sum(a[(i, j)] for j in neighb[i]) <= 1)  # Bornitude
    @constraint(
        model,
        exit_region[r in keys(c_data.regions_to_nodes)],
        sum(
            a[(i, j)] for i in c_data.regions_to_nodes[r]
            for j in setdiff(neighb[i], c_data.regions_to_nodes[r])
        ) >= 1,
    )  # Sortie de chacune des regions
    @constraint(model, min_visits, sum(a) >= c_data.Amin - 1)  # Visite au moins Amin aerodromes
    elem_constraints = Vector{Vector{Any}}(undef, 0)  # Assure que le chemin est elementaire
    @objective(model, Min, sum(c_data.D[i, j] * a[(i, j)] for (i, j) in arcs))
    unset_binary.(a)

    violated_constraints = true
    tilde_a = Dict{Tuple{Int64,Int64},Float64}()
    while violated_constraints
        optimize!(model)
        println(termination_status(model))
        @assert termination_status(model) == OPTIMAL
        @assert primal_status(model) == FEASIBLE_POINT
        for (i, j) in arcs
            tilde_a[(i, j)] = value(a[(i, j)])
        end
        v_value, v_set = separation_pb(pb_data=pb_data, tilde_a=tilde_a)
        println("violation value = ", v_value)
        v_arcs = Set{Tuple{Int64,Int64}}((i, j) for (i, j) in arcs if i in v_set && j in v_set)
        if v_value < epsilon
            if all(is_binary.(a))
                violated_constraints = false
            else
                set_binary.(a)
            end
        else
            new_constraints = @constraint(
                model,
                [k in v_set],
                sum(a[(i, j)] for (i, j) in v_arcs) <= sum(a[(i, j)] for i in v_set for j in neighb[i] if i != k),
            )
            push!(elem_constraints, new_constraints)
        end
    end
    chosen_arcs = Set{Tuple{Int64,Int64}}((i, j) for (i, j) in arcs if value(a[(i, j)]) > 0.5)
    path_length = sum(c_data.D[i, j] for (i, j) in chosen_arcs)
    path = Path(makepath(chosen_arcs), path_length)
end




#### Fonctions de plot





function edge_matrix(m)
    distance_matrix = Dict{Tuple{Int,Int},Float64}()
    for i in 1:size(m, 1)
        for j in 1:size(m, 2)
            distance_matrix[(i, j)] = m[i, j]
        end
    end
    return distance_matrix
end


function distance_matrix_to_edges_and_weights(distance_matrix)
    edges = Tuple{Int,Int}[]  # Initialize an empty array for edges
    weights = Float64[]       # Initialize an empty array for edge weights

    for i in 1:size(distance_matrix, 1)
        for j in 1:size(distance_matrix, 2)
            if distance_matrix[i, j] != 0.0
                push!(edges, (i, j))
                push!(weights, distance_matrix[i, j])
            end
        end
    end

    return edges, weights
end

function create_coords(matrix)
    x = matrix[:, 1]
    y = matrix[:, 2]
    return x, y
end

function make_path(nodes)
    path = []
    for i in 1:length(nodes)-1
        push!(path, (nodes[i], nodes[i+1]))
    end
    return path
end

function nodes_in_order(path)
    node_list = Int[]

    for edge in path
        n1, n2 = edge
        push!(node_list, n1)
        push!(node_list, n2)
    end

    unique_nodes = Int[]
    [push!(unique_nodes, node) for node in node_list if node ∉ unique_nodes]

    return unique_nodes
end


function max_node(graph)
    max_node = -Inf

    for edge in graph
        max_node = max(max_node, maximum(edge))
    end

    return max_node
end


function convert_rgb(seq)
    conv = [convert(RGB{FixedPointNumbers.Normed{UInt8,8}}, rgb) for rgb in seq]
    return conv
end

function color_node(path, nb_nodes)
    start_color = RGB(1.0, 0.0, 0.0)  # Red
    end_color = RGB(0.5, 0.0, 0.5)    # Purple

    length_path = length(path)
    node_color = []
    col_index = 1

    # Create a dictionary to map nodes to their positions in the path
    node_positions = Dict(node => pos for (pos, node) in enumerate(path))

    for i in 1:nb_nodes
        if i in keys(node_positions)
            # Determine the color based on the position in the path
            position = node_positions[i]
            color = start_color + (position - 1) / (length_path - 1) * (end_color - start_color)
            push!(node_color, color)
        else
            push!(node_color, RGB(0.5, 0.5, 0.5))
        end
    end

    return node_color
end


function plot_path(edges1, edges2, weight2, coords, name)

    vec_xNode, vec_yNode = create_coords(coords)

    gr = SimpleDiGraph()
    Graphs.add_vertices!(gr, length(vec_xNode))
    weights_dict = Dict(zip(edges2, weight2))

    for edge in edges2
        Graphs.add_edge!(gr, edge)
    end

    mgr = MetaDiGraph(gr)

    # Add the attribute of nodes
    for i in 1:length(vec_xNode)
        set_props!(mgr, i, Dict(
            Symbol("vec_xNode") => vec_xNode[i],
            Symbol("vec_yNode") => vec_yNode[i],
        ))
    end

    color_nodes = convert_rgb(color_node(nodes_in_order(edges2), max_node(edges1)))
    println(color_nodes)


    gp = graphplot(
        gr,
        x=vec_xNode,
        y=vec_yNode,
        names=1:Graphs.nv(gr),
        fontsize=8,
        nodeshape=:circle,
        markersize=10,
        markerstrokewidth=1,
        # edges
        edgelabel=weights_dict,  # Use the edge_labels array for labels
        edgelabelfontsize=1,
        edgelabelposition="above",
        linewidth=2,
        markercolor=color_nodes,
        curvature_scalar=0.2,
    )
    savefig(gp, name)
end






#### Solve + plot



instance = joinpath("data", "instance_20_1.txt")
n, d, f, Amin, Nr, R, regions, coord, D = readInstance(instance)
pb_data = create_pb_data(instance)
c_data = pb_data.c_data

println("nb sommets = ", c_data.n)
println("depart = ", c_data.d)
println("fin = ", c_data.f)
println("min nb aerodromes a visiter = ", c_data.Amin)
println("nb regions = ", c_data.Nr)

path = expo_solve(pb_data)


#plot
x, y = create_coords(coord)
graph1, weight1 = distance_matrix_to_edges_and_weights(D)


node = path.nodes

node_vect = []
for i in node
    push!(node_vect, i)
end

graph2 = make_path(node_vect)

weight2 = []
for edge in graph2
    push!(weight2, edge_matrix(D)[edge])
end

plot_path(graph1, graph2, weight2, coord, "inst_20_1_opti.png")