module GraphingPlanes

using DataFrames
using GeoMakie, CairoMakie
using GraphMakie, GraphMakie.Graphs
using NaturalEarth

using Graphs, GraphPlot
using SimpleWeightedGraphs
using Colors


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
  positions = Dict(row.row_id => [row.LONGITUDE, row.LATITUDE] for row in eachrow(US_coords))
  positions_vector = Point2f.(collect(values(positions)))

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



function us_graph(US_data::DataFrame, US_coords::DataFrame)
  positions = Dict(row.row_id => [row.LONGITUDE, row.LATITUDE] for row in eachrow(US_coords))
  positions_vector = Point2f.(collect(values(positions)))

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
    add_edge!(g, origin_id, dest_id)
  end

  println(size(g))

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
    layout = _ -> positions_vector, node_size = 5,
    edge_color = (:blue, 0.01)
    # TODO: Define edge color by degree or some other factor
    # TODO: Figure out transparency for a range of colors
    # edge_color = cgrad(:plasma, 0.05)[LinRange(0, 1, ne(g))]
  )

  return fig
end





end