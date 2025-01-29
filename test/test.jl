using OSMToolset
import OpenStreetMapX
import OpenStreetMapX: ENU
using DataFrames
using CSV
using Test
include("../src/OSMXGraph.jl")

dir_in = "./example/data"

road_types = ["motorway", "trunk", "primary", "secondary", 
            "tertiary", "residential", "service", "living_street", 
            "motorway_link", "trunk_link", "primary_link", "secondary_link", 
            "tertiary_link"] 

osm_file = "Ochota.osm"


parsed = OpenStreetMapX.parseOSM(string(dir_in,"/",osm_file))
ways = parsed.ways
filtered_ways = OSMXGraph.filter_ways(ways,road_types)

tst = [(haskey(way.tags, "highway") && (way.tags["highway"] in road_types)) for way in ways]
tst_1 = findall(x -> x == 1, tst)

@test length(filtered_ways) == length(tst_1)

ways_intersections, intersections, road_tags, nodes = OSMXGraph.find_intersections(filtered_ways, parsed)

@test length(ways_intersections) == 5995
@test length(intersections) == 8108
@test length(road_tags) == 5995
@test length(nodes) == 8108

edges = OSMXGraph.ways_to_edges(ways_intersections,road_tags,parsed,nodes)

@test length(edges) == 16041

df = OSMXGraph.edges_to_df(edges)

@test size(df,1) == 16041
@test size(df,2) == 9

sparse_index = OSMXGraph.create_sparse_index(df.from_id,df.to_id,df.id)

@test size(sparse_index,1) == 8108
@test size(sparse_index,2) == 8108

df, sparse_index, road_index, node_ids = OSMXGraph.create_road_graph("Ochota.osm", road_types,"Ochota_graph.csv","Ochota_nodes.json",dir_in=dir_in)
POI_df = OSMToolset.find_poi(string(dir_in,"/","Ochota.osm"))

@test size(POI_df,2) == 10

POI_xs = POI_df.lat
POI_ys = POI_df.lon
poi_with_nearest_points = OSMXGraph.add_nearest_road_point(POI_df,POI_xs, POI_ys, road_index, node_ids)

@test size(poi_with_nearest_points,2) == 11