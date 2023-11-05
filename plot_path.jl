using Graphs, GraphRecipes, Plots, Colors, GraphPlot, FixedPointNumbers,MetaGraphs



distance_matrix = Dict{Tuple{Int, Int}, Float64}()
for i in 1:size(m, 1)
    for j in 1:size(m, 2)
        distance_matrix[(i, j)] = m[i, j]
    end
end
return distance_matrix
end


function distance_matrix_to_edges_and_weights(distance_matrix)
edges = Tuple{Int, Int}[]  # Initialize an empty array for edges
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


function euclidean_distance_matrix(x, y)
n = length(x)
m = length(y)

distance_matrix = zeros(n, m)

for i in 1:n
    for j in 1:m
        distance_matrix[i, j] = round(sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2))
    end
end

return distance_matrix
end



function create_coords(matrix)
x = matrix[:, 1]
y = matrix[:, 2]
return x, y
end

function make_path(nodes)
path = []
for i in 1:length(nodes) - 1
    push!(path, (nodes[i], nodes[i + 1]))
end
return path
end

function nodes(path)
node_set = Set{Int}()

for edge in path
    n1, n2 = edge
    push!(node_set, n1)
    push!(node_set, n2)
end

return collect(node_set)
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


function make_graph(edges,weights)
    
g = Graphs.DiGraph()
n_edges = length(edges)
n_max = max_node(edges)
edgelabels = Dict{Tuple{Int, Int}, Float64}() 
    
for j in 1:n_max
    Graphs.add_vertices!(g, 1)
end

for k in 1:n_edges
    ps, pe = edges[k]
    Graphs.add_edge!(g, ps, pe)
    edgelabels[(ps, pe)] = weights[k]  
end

return g,edgelabels
end
    
function directed_adjacency_matrix(edges, num_nodes)

adj_matrix = zeros(Int, num_nodes, num_nodes)

for (src, dest) in edges
    adj_matrix[src, dest] = 1
end

return adj_matrix
end

function convert_rgb(seq)
conv = [convert(RGB{FixedPointNumbers.Normed{UInt8, 8}}, rgb) for rgb in seq]
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

function plot_path4(edges1,edges2,weight2, coords,name)

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


gp= graphplot(
    gr,
    x = vec_xNode,
    y = vec_yNode,
    names = 1:Graphs.nv(gr),
    fontsize = 8,
    nodeshape = :circle,
    markersize = 1,
    markerstrokewidth = 1,
    # edges
    edgelabel = weights_dict,  # Use the edge_labels array for labels
    edgelabelfontsize = 1,
    edgelabelposition = "above",
    linewidth = 2,
    markercolor = color_nodes,
    curvature_scalar = 0.2,
)

savefig(gp, name)
end


function plot_path(edges1,edges2,weight2, coords,name)

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


graphplot(
    gr,
    x = vec_xNode,
    y = vec_yNode,
    names = 1:Graphs.nv(gr),
    fontsize = 8,
    nodeshape = :circle,
    markersize = 1,
    markerstrokewidth = 1,
    # edges
    edgelabel = weights_dict,  # Use the edge_labels array for labels
    edgelabelfontsize = 1,
    edgelabelposition = "above",
    linewidth = 2,
    markercolor = color_nodes,
    curvature_scalar = 0.2,
)
end


#solve the prblm
instance = joinpath("data", "instance_40_1.txt")
n, d, f, Amin, Nr, R, regions, coord, D  = readInstance(instance)
pb_data = create_pb_data(instance)
c_data = pb_data.c_data

println("nb sommets = ", c_data.n)
println("depart = ", c_data.d)
println("fin = ", c_data.f)
println("min nb aerodromes a visiter = ", c_data.Amin)
println("nb regions = ", c_data.Nr)

poly_solve(pb_data)


#plot
x,y = create_coords(coord)
matrix = euclidean_distance_matrix(x, y)
graph1,weight1 = distance_matrix_to_edges_and_weights(matrix)


node = paths.nodes

node_vect = []
for i in node
    push!(node_vect,i)
end

graph2 = make_path(node_vect)

weight2 = []
for edge in graph2
    push!(weight2,edge_matrix(matrix)[edge])
end

plot_path(graph1,graph2,weight2, coord, "inst_40_1_opti.png")