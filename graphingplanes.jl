module GraphingPlanes

using DataFrames
using GeoMakie, CairoMakie
using GraphMakie, GraphMakie.Graphs
using NaturalEarth

using Graphs, GraphPlot
using SimpleWeightedGraphs
using Colors
import Statistics


# Remove invalid edges (out-of-range endpoints) and self-loops from graphs so plotting positions align
function sanitize_graph(g::SimpleWeightedDiGraph)
  n = nv(g)
  g2 = SimpleWeightedDiGraph(n)
  w = weights(g)
  removed = 0
  for e in edges(g)
    if 1 <= e.src <= n && 1 <= e.dst <= n && e.src != e.dst
      add_edge!(g2, e.src, e.dst, w[e.src, e.dst])
    else
      removed += 1
    end
  end
  if removed > 0
    @info "sanitize_graph: removed $removed invalid/self-loop edges from SimpleWeightedDiGraph"
  end
  return g2
end

function sanitize_graph(g::SimpleGraph)
  n = nv(g)
  g2 = SimpleGraph(n)
  removed = 0
  for e in edges(g)
    if 1 <= e.src <= n && 1 <= e.dst <= n && e.src != e.dst
      add_edge!(g2, e.src, e.dst)
    else
      removed += 1
    end
  end
  if removed > 0
    @info "sanitize_graph: removed $removed invalid/self-loop edges from SimpleGraph"
  end
  return g2
end


### Helper Functions ###
function standardize(degrees, max_val)
  dmin = minimum(degrees)
  dmax = maximum(degrees)
  if dmax == dmin
    # all degrees equal -> use the midpoint of the target range
    scaled = fill((5.0 + max_val) / 2, length(degrees))
  else
    scaled = 5.0 .+ (float.(degrees) .- dmin) .* (max_val - 5.0) ./ (dmax - dmin)
  end
  return convert.(Int, round.(scaled))
end

##### Networking Functions #####

#Helper functions to compute in and out degrees
function weighted_outdegree(g::SimpleWeightedDiGraph)
    w = weights(g)
    [sum(w[i, neighbors(g, i)]) for i in 1:nv(g)]
end

function weighted_indegree(g::SimpleWeightedDiGraph)
    w = weights(g)
    [sum(w[inneighbors(g, i), i]) for i in 1:nv(g)]
end



function FindHubs(g::SimpleWeightedDiGraph, id_to_airport::Dict)
    #compute degrees
    outdeg = weighted_outdegree(g)
    indeg = weighted_indegree(g)
    #find total connections
    totaldeg = outdeg .+ indeg

    #Sort the nodes by total degree in descending order
    sorted_nodes = partialsortperm(totaldeg, 1:5, rev=true)
    top_airports = [(id_to_airport[i], totaldeg[i]) for i in sorted_nodes]

    return top_airports
end

function FindBridges(g::SimpleWeightedDiGraph, id_to_airport::Dict; top_n=5)
    # Compute betweenness centrality (unweighted)
    bc = betweenness_centrality(g)

    # Find top nodes
    top_idx = partialsortperm(bc, 1:top_n, rev=true)
    top_airports = [(id_to_airport[i], bc[i]) for i in top_idx]

    return top_airports
end


function FindBridges(g::SimpleGraph, id_to_airport::DataFrame; top_n=5)
    # Compute betweenness centrality (unweighted)
    bc = betweenness_centrality(g)

    # Find top nodes
    top_idx = partialsortperm(bc, 1:top_n, rev=true)
    top_airports = [(id_to_airport[i], bc[i]) for i in top_idx]

    return top_airports
end



function BuildGraph(df::DataFrame)
#GroupBy origin airport code and destination airport code
# this gives unique pairs and the number of passengers will be summed to be used as the edge's weight

    grouped = groupby(df, [:OriginAirportAlpha, :DestinationAirportAlpha])
    edges = combine(grouped, :PassengersTransported => sum => :weight)
    #Ignore Freight planes
    edges = filter(:weight => w -> w > 0, edges)

    # Removes Freight Planes from Df
    filter!(:DestinationAirportAlpha => A -> A in edges.DestinationAirportAlpha, df)
    filter!(:OriginAirportAlpha => A -> A in edges.OriginAirportAlpha, df)

    #Define a vector of unique nodes
    nodes = union(df[!, :OriginAirportAlpha], df[!, :DestinationAirportAlpha])

    #map nodes to codes (ID's)
    airport_to_id = Dict(nodes[i] => i for i in eachindex(nodes))
    id_to_airport = Dict(v => k for (k, v) in airport_to_id)

    #empty graph
    g = SimpleWeightedDiGraph(length(nodes))

    #Add weighted edges
    for row in eachrow(edges)
        origin_id = airport_to_id[row.OriginAirportAlpha]
        dest_id = airport_to_id[row.DestinationAirportAlpha]
        weight = row.weight
        add_edge!(g, origin_id, dest_id, weight)
    end

    airport_labels = [id_to_airport[i] for i in 1:nv(g)]
    gplot(g, nodelabel=airport_labels)
    return g, id_to_airport
end

function BuildUSDiGraph(df::DataFrame, US_coords::DataFrame)
    # Filter to flights with actual passengers
    grouped = groupby(df, [:Origin, :Dest])
    edges = combine(grouped, :PassengersTransported => sum => :weight)
    edges = filter(:weight => w -> w > 0, edges)

    # Build a clean mapping: airport name => node ID (1-based index in filtered_coords)
    airport_to_id = Dict(row.DisplayName => i for (i, row) in enumerate(eachrow(US_coords)))
    id_to_airport = Dict(v => k for (k, v) in airport_to_id)

    # Create the graph with as many nodes as filtered_coords
    g = SimpleWeightedDiGraph(nrow(US_coords))

    #Add weighted edges
    for row in eachrow(edges)
        origin_id = airport_to_id[row.Origin]
        dest_id = airport_to_id[row.Dest]
        weight = row.weight
        add_edge!(g, origin_id, dest_id, weight)
    end


    return g, id_to_airport
end




function BuildUSGraph(df::DataFrame, US_coords::DataFrame)
  g = SimpleGraph(nrow(US_coords))

  for row in eachrow(df)
    if row.Origin in US_coords.DisplayName
      origin_id = US_coords[US_coords.DisplayName .== row.Origin, :row_id][1]
    else
      continue
    end
    if row.Dest in US_coords.DisplayName
      dest_id = US_coords[US_coords.DisplayName .== row.Dest, :row_id][1]
    else
      continue
    end
    add_edge!(g, origin_id, dest_id)
  end

  return g
end



function PlotGraph(g::SimpleWeightedDiGraph, id_to_airport::Dict; title="Season Network")
    #Invert the dictionary to get airports to node ID
    airport_to_id = Dict(v => k for (k, v) in id_to_airport)

    #use the top 5 airports
    top_airports = FindHubs(g, id_to_airport)
    top_ids = [airport_to_id[a[1]] for a in top_airports]

    #Build a labels list
    airport_labels = [id_to_airport[i] for i in 1:nv(g)]

    #Define colors, highlight the hubs
    colors = fill(colorant"lightblue", nv(g))
    for i in top_ids
        colors[i] = colorant"orange"
    end

    #Define node sizes proportional to total degree
    w = weights(g)
    sizes = [sum(w[i, :]) + sum(w[:, i]) for i in 1:nv(g)]
    sizes = sizes ./ maximum(sizes) .* 15  # scale for visibility


    # Plot the graph
    gplot(g;
        nodelabel=airport_labels,
        nodefillc=colors,
        nodesize=sizes,
        arrowlengthfrac=0.05,
        title=title
    )
end






##### Geo Mapping #####

admin_1_df = DataFrame(naturalearth("admin_1_states_provinces", 110))
filter!(:gu_a3 => ==("USA"), admin_1_df)

admin_1_df.geometry = GeoMakie.to_multipoly.(admin_1_df.geometry)
num_countries = size(admin_1_df, 1)


function us_pointmap(US_data::DataFrame, US_coords::DataFrame)
  # Build positions vector in the same order as US_coords row indices so node i maps to US_coords row i
  positions_vector = [Point2f(US_coords[i, :LONGITUDE], US_coords[i, :LATITUDE]) for i in 1:nrow(US_coords)]

  g = SimpleGraph(nrow(US_coords))

  for row in eachrow(US_data)
    if row.Origin in US_coords.DisplayName
      origin_id = US_coords[US_coords.DisplayName .== row.Origin, :row_id][1]
    else
      continue
    end
    if row.Dest in US_coords.DisplayName
      dest_id = US_coords[US_coords.DisplayName .== row.Dest, :row_id][1]
    else
      continue
    end
    # add_edge!(g, origin_id, dest_id)
  end


  fig = Figure(size = (1200, 800), fontsize = 22)
  ga = GeoAxis(
      fig[1, 1],
      source = "+proj=longlat +datum=WGS84",
      dest = "+proj=lcc +lon_0=-100 +lat_1=33 +lat_2=45",
      title = "US Airline Locations",
  )

  poly!(
    ga, admin_1_df.geometry;
    color = 1:num_countries, colormap = (:viridis, 0.25),
    strokecolor = :black, strokewidth = 1
  )

  graphplot!(
    ga, g;
    layout = _ -> positions_vector, node_size = 10
  )

  return fig
end




function us_graph(g::SimpleWeightedDiGraph, US_coords::DataFrame, title::String="US Airline Network")
  # Positions must be in node index order (1..nrow(US_coords)) so graph nodes align with coordinates
  positions_vector = [Point2f(US_coords[i, :LONGITUDE], US_coords[i, :LATITUDE]) for i in 1:nrow(US_coords)]

  # sanitize graph to remove self-loops or edges with invalid node indices
  g = sanitize_graph(g)

  if length(positions_vector) != nv(g)
      @warn "positions vector length ($(length(positions_vector))) != number of nodes ($(nv(g)))."
  end

  fig = Figure(size = (1200, 800), fontsize = 22)
  ga = GeoAxis(
      fig[1, 1],
      source = "+proj=longlat +datum=WGS84",
      dest = "+proj=lcc +lon_0=-100 +lat_1=33 +lat_2=45",
      title = title,
  )

  poly!(
    ga, admin_1_df.geometry;
    color = 1:num_countries, colormap = (:viridis, 0.25),
    strokecolor = :black, strokewidth = 1
  )

  scaler = standardize(degree(g), 25)
  top1 = Statistics.quantile(scaler, 0.99)
  top25 = Statistics.quantile(scaler, 0.75)
  bottom25 = Statistics.quantile(scaler, 0.25)

  # Map values to categories: 1 (bottom 25%), 2 (middle), 3 (top 25%)
  membership = map(x -> x >= top1 ? 4 : x >= top25 ? 3 : x <= bottom25 ? 1 : 2, scaler)

  nodecolor = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=4)

  # membership color
  nodefillc = nodecolor[membership]

  # Edge Colors based on weights
  # 1. Get edge weights
  edges_list = collect(edges(g))
  weights = [g.weights[e.src, e.dst] for e in edges_list]

  # 2. Quantiles for grouping (25%, 75%, etc.)
  sorted_weights = sort(weights)
  bottom25 = Statistics.quantile(sorted_weights, 0.25)
  top25    = Statistics.quantile(sorted_weights, 0.75)
  top1     = Statistics.quantile(sorted_weights, 0.99)

  # 3. Bucket edges by weight
  edge_membership = map(x -> x >= top1 ? 4 : x >= top25 ? 3 : x <= bottom25 ? 1 : 2, weights)

  # 4. Create color gradient (same as node color)
  base_colors = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=4)
  base_colors_rgb = convert.(RGB, base_colors)  # HSL → RGB

  # # 5. Normalize weights for alpha (0.1 for lightest, up to 1.0 for heaviest)
  # min_w, max_w = minimum(weights), maximum(weights)
  # alpha_vals = 0.1 .+ 0.9 .* ((weights .- min_w) ./ (max_w - min_w + eps()))  # scaled to [0.1, 1.0]

  # 5. Assign alpha: 1.0 if in top 1%, else 0.0 (fully transparent)
  alpha_vals = map(w -> w >= top1 ? 1.0 : 0.0, weights)

  # 6. Combine color and alpha
  edge_colors = RGBA.(base_colors_rgb[edge_membership], alpha_vals)


  graphplot!(
    ga, g;
    layout = _ -> positions_vector, node_size = scaler,
    # nodelabel = string.(US_coords.DisplayName), nodelabelsize = 10,
    node_color = nodefillc,
    edge_color = edge_colors,
    markersize=0
    # edge_color = cgrad(:plasma, 0.05)[LinRange(0, 1, ne(g))]
  )

  return fig
end




function us_graph(US_data::DataFrame, US_coords::DataFrame, title::String="US Airline Network")
  # Positions must be in node index order (1..nrow(US_coords)) so graph nodes align with coordinates
  positions_vector = [Point2f(US_coords[i, :LONGITUDE], US_coords[i, :LATITUDE]) for i in 1:nrow(US_coords)]

  g = BuildUSGraph(US_data, US_coords)

  # sanitize graph to remove self-loops or edges with invalid node indices
  g = sanitize_graph(g)

  if length(positions_vector) != nv(g)
      @warn "positions vector length ($(length(positions_vector))) != number of nodes ($(nv(g)))."
  end

  fig = Figure(size = (1200, 800), fontsize = 22)
  ga = GeoAxis(
      fig[1, 1],
      source = "+proj=longlat +datum=WGS84",
      dest = "+proj=lcc +lon_0=-100 +lat_1=33 +lat_2=45",
      title = title,
  )

  poly!(
    ga, admin_1_df.geometry;
    color = 1:num_countries, colormap = (:viridis, 0.25),
    strokecolor = :black, strokewidth = 1
  )

  scaler = standardize(degree(g), 25)
  top1 = Statistics.quantile(scaler, 0.99)
  top25 = Statistics.quantile(scaler, 0.75)
  bottom25 = Statistics.quantile(scaler, 0.25)

  # Map values to categories: 1 (bottom 25%), 2 (middle), 3 (top 25%)
  membership = map(x -> x >= top1 ? 4 : x >= top25 ? 3 : x <= bottom25 ? 1 : 2, scaler)

  nodecolor = range(HSL(colorant"red"), stop=HSL(colorant"green"), length=4)

  # membership color
  nodefillc = nodecolor[membership]

  graphplot!(
    ga, g;
    layout = _ -> positions_vector, node_size = scaler,
    nodelabel = string.(US_coords.DisplayName), nodelabelsize = 10,
    node_color = nodefillc,
    edge_color = (:blue, 0.01)
    # edge_color = cgrad(:plasma, 0.05)[LinRange(0, 1, ne(g))]
  )

  return fig
end

function us_graph(US_data::DataFrame, US_coords::DataFrame, communities::Vector{BitSet}, title::String="US Airline Network")
  # Positions must be in node index order (1..nrow(US_coords)) so graph nodes align with coordinates
  positions_vector = [Point2f(US_coords[i, :LONGITUDE], US_coords[i, :LATITUDE]) for i in 1:nrow(US_coords)]

  g = BuildUSGraph(US_data, US_coords)

  # sanitize graph to remove self-loops or edges with invalid node indices
  g = sanitize_graph(g)

  if length(positions_vector) != nv(g)
      @warn "positions vector length ($(length(positions_vector))) != number of nodes ($(nv(g)))."
  end

  fig = Figure(size = (1200, 800), fontsize = 22)
  ga = GeoAxis(
      fig[1, 1],
      source = "+proj=longlat +datum=WGS84",
      dest = "+proj=lcc +lon_0=-100 +lat_1=33 +lat_2=45",
      title = "US Airline Locations",
  )

  poly!(
    ga, admin_1_df.geometry;
    color = 1:num_countries, colormap = (:viridis, 0.25),
    strokecolor = :black, strokewidth = 1
  )


  membership = zeros(Int, nv(g))

  for (comm_id, bitset) in enumerate(communities)
      for node in bitset
          membership[node] = comm_id
      end
  end

  membership .+= 1

  nodecolor = range(HSV(0,1,1), stop=HSV(-360,1,1), length=maximum(membership))

  # membership color
  nodefillc = nodecolor[membership]


  scaler = standardize(degree(g), 25)
  top1 = Statistics.quantile(scaler, 0.99)
  top25 = Statistics.quantile(scaler, 0.75)
  bottom25 = Statistics.quantile(scaler, 0.25)

  graphplot!(
    ga, g;
    layout = _ -> positions_vector, node_size = scaler,
    # nodelabel = string.(US_coords.DisplayName), nodelabelsize = 10,
    node_color = nodefillc,
    edge_color = (:blue, 0.01)
    # edge_color = cgrad(:plasma, 0.05)[LinRange(0, 1, ne(g))]
  )

  return fig
end


### Diagnostics helpers
"""
Find nodes in `US_coords` that are geographically isolated: their nearest neighbor distance
is greater than `radius_deg` (in degrees). Returns a DataFrame with row_id, DisplayName, LONGITUDE, LATITUDE, nearest_dist.
"""
function find_isolated_nodes(US_coords::DataFrame; radius_deg::Float64=5.0)
  n = nrow(US_coords)
  if n < 2
    return DataFrame()
  end
  lons = Float64.(US_coords.LONGITUDE)
  lats = Float64.(US_coords.LATITUDE)

  nearest_dists = fill(Inf, n)
  for i in 1:n
    for j in 1:n
      if i == j
        continue
      end
      dlon = lons[i] - lons[j]
      dlat = lats[i] - lats[j]
      dist = sqrt(dlon^2 + dlat^2)
      if dist < nearest_dists[i]
        nearest_dists[i] = dist
      end
    end
  end

  rows = Vector{NamedTuple}()
  for i in 1:n
    if nearest_dists[i] > radius_deg
      push!(rows, (row_id=US_coords.row_id[i], DisplayName=US_coords.DisplayName[i], LONGITUDE=US_coords.LONGITUDE[i], LATITUDE=US_coords.LATITUDE[i], nearest_dist=nearest_dists[i]))
    end
  end
  return DataFrame(rows)
end

"""
List edges in a graph `g` that are incident to node `node_id` (1-based). Returns a vector of (src,dst) pairs.
"""
function edges_for_node(g, node_id::Integer)
  return [(e.src, e.dst) for e in edges(g) if e.src == node_id || e.dst == node_id]
end

"""
Run quick diagnostics: print isolated nodes in `US_coords` (default radius 5°) and for each such node, print incident edges in `g`.
"""
function diagnose_misplaced_nodes(g, US_coords::DataFrame; radius_deg::Float64=5.0)
  iso = find_isolated_nodes(US_coords; radius_deg=radius_deg)
  if nrow(iso) == 0
    println("No isolated nodes found with radius_deg=$(radius_deg). Try lowering radius.")
    return iso
  end

  println("Found $(nrow(iso)) isolated nodes (nearest neighbor > $(radius_deg) degrees):")
  show(stdout, iso)

  for r in eachrow(iso)
    node = r.row_id
    println("\nEdges incident to node $(node) ($(r.DisplayName)):")
    println(edges_for_node(g, node))
  end
  return iso
end

end