## Tools for transforming OSM data into a graph structure

The 'OSMXGraph' module was created under contract no. 7474/2024. The tools enable the transformation of OSM data into a set of data frames containing:
1. A graph of the road network (detailed in the "Road Network Graph" section),
2. Metadata for the graph's edges (detailed in the "Dataframe for Graph Metadata" section),
3. Metadata for the graph's nodes (detailed in the "Dataframe for Graph Metadata" section),s
4. Metadata for POIs, including the point's location, class and type of POI, and the ID of the nearest road nodes (detailed in the "Dataframe for POI" section).

## Road Network Graph

The road network graph is stored as a sparse adjacency matrix, where each pair of node IDs is mapped to a corresponding road ID. The node and way IDs can then be mapped to a dataframe containing metadata about the nodes and roads.

## Dataframe for Graph Metadata

The dataframe contains metadata for both road nodes and roadways. Each row of the structure represents a one-way road. The dataframe has the following columns:

- **id**: Road's ID.
- **from_id**: ID of the node that is the starting point of the road. The ID corresponds to the node ID in the road network graph and ID in the POI dataframe.
- **to_id**: ID of the node that is the endpoint of the road. The ID corresponds to the node ID in the road network graph and ID in the POI dataframe.  
- **from**: ID of the node that is the starting point of the road. The ID corresponds to the node ID in the OSM file.
- **to**: ID of the node that is the endpoint of the road. The ID corresponds to the node ID in the OSM file.
- **from_LLA**: Location of the starting point of the road in geodetic metric.
- **to_LLA**: Location of the endpoint of the road in geodetic coordinates.
- **way**: Road's ID in the OSM file.
- **type**: Road's type.

Information about edges and nodes is also available as two separate data frames.
These data frames contain the same information as the data frame for graph metadata but are split into two independent structures. The 'way' column in the metadata data frame enables mapping each road to its corresponding metadata dictionary

## Dataframe for POI

The dataframe contains locations of points of interest, their metadata, and the location of the nearest road node. The structure is based on the OSMToolset POI dataframe. The module extends this structure by adding the ID of the nearest road's node.

## Overview of the Graph Building Process


```julia
using OSMXGraph
using OSMToolset
using OpenStreetMapX

dir_in = "../data"
road_types = ["motorway", "trunk", "primary", "secondary", 
            "tertiary", "residential", "service", "living_street", 
            "motorway_link", "trunk_link", "primary_link", "secondary_link", 
            "tertiary_link"] 
osm_file = "Warszawa.osm"
graph_file_name = "Warszawa_graph.csv"
node_file_name = "Warszawa_nodes.json"
dir_in=dir_in
```


### 1. Filtering ways
The First step in creating graph is filtering the vector of 'Way' objects. The 'Way' 
structure in OSM files is a vector that represents more types of objects than just roads.
The 'filter ways' function filters the vector of 'Way' objects to retain only roads with 
hierarchy types specified by user.  


```julia
parsed = OpenStreetMapX.parseOSM(string(dir_in,"/",osm_file))
ways = parsed.ways
filtered_ways = OSMXGraph.filter_ways(ways,road_types)
```


    125255-element Vector{Way}:
     Way(4307329, [2448759046, 7093785352, 2452307268, 1439574696], Dict("name:etymology:wikidata" => "Q5441838", "surface" => "asphalt", "name" => "Rondo Feliksa Stamma", "sidewalk:right" => "separate", "wikidata" => "Q113528575", "lit" => "yes", "highway" => "tertiary", "cycleway:both" => "no", "junction" => "roundabout", "sidewalk:left" => "no"…))
     Way(4307330, [26063923, 26063924, 9900672046, 7093724913, 26063925, 4770319996, 7554494461, 9442368604, 11303130233, 11303232206, 9220363525, 9982770288], Dict("lit" => "yes", "name" => "Bokserska", "source:maxspeed" => "PL:urban", "highway" => "tertiary", "lanes" => "2", "maxspeed" => "50", "surface" => "asphalt", "wikidata" => "Q106807412"))
     Way(4308966, [3387797238, 9252948452], Dict("cycleway:left" => "no", "cycleway:right" => "separate", "lit" => "yes", "name" => "Postępu", "highway" => "tertiary", "lanes" => "2", "sidewalk" => "separate", "maxspeed" => "50", "surface" => "asphalt", "wikidata" => "Q63188829"…))
     Way(4308978, [3576562134, 8961134644, 9936189679, 4972766762, 26083936, 3705172931, 3576562138, 9035856380, 427569388, 9035856661, 427569361, 427569362, 260821813, 427569360, 9035856642, 260821818, 9035856643, 9035856663, 4972766763, 260821819], Dict("lit" => "yes", "highway" => "residential", "surface" => "asphalt", "noname" => "yes"))
     Way(4308979, [26083933, 4439989380, 3024770845, 8677742938, 9281835433, 9934684736, 10692999753, 1147905901, 1147905852, 8975524567, 26083935], Dict("highway" => "residential", "sidewalk" => "separate", "surface" => "asphalt", "noname" => "yes"))
     Way(4308980, [26083939, 427569365, 3024770846, 9281835907, 26083947], Dict("lit" => "yes", "name" => "Jana Zaorskiego", "sidewalk:left" => "no", "oneway" => "yes", "sidewalk:right" => "separate", "highway" => "residential", "lanes" => "1", "surface" => "asphalt"))


### 2. Finding nodes and edges
The step is performed based on filtered roadway vector. The function 'find intersaction' 
iterates through all roads and their nodes. A road point is considered a graph node 
if it is a starting or ending point of a road or if the it is an intersecion of 
roads (i.e, the function has already encountered the point during the iteration). As a result, 
each roadway retains only those points that are starting/ending points or 
intersections. Edges are created in the 'ways\_to\_edges' function by pairing the nearest 
nodes. Edges represent one-way roads, so if a road is bidirectional, the roadway is 
split into one-directional ways. 


```julia
ways_intersections, intersections, road_tags, nodes = OSMXGraph.find_intersections(filtered_ways, parsed)
edges = OSMXGraph.ways_to_edges(ways_intersections,road_tags,parsed,nodes)
```


    362340-element Vector{Main.OSMXGraph.Edge}:
     Main.OSMXGraph.Edge(1, 63583, 2053, 4978125632, 1181144372, LLA(52.2566118, 21.0331274, 0.0), LLA(52.2566598, 21.0331933, 0.0), 451965551, "tertiary")
     Main.OSMXGraph.Edge(2, 34902, 34901, 1949648310, 9265228166, LLA(52.2159197, 20.979194, 0.0), LLA(52.2154663, 20.9783547, 0.0), 184475036, "service")
     Main.OSMXGraph.Edge(3, 34901, 34902, 9265228166, 1949648310, LLA(52.2154663, 20.9783547, 0.0), LLA(52.2159197, 20.979194, 0.0), 184475036, "service")
     Main.OSMXGraph.Edge(4, 89783, 10564, 8213623939, 1495256494, LLA(52.2757071, 21.0658142, 0.0), LLA(52.2767636, 21.0666679, 0.0), 29448268, "service")
     Main.OSMXGraph.Edge(5, 10564, 89783, 1495256494, 8213623939, LLA(52.2767636, 21.0666679, 0.0), LLA(52.2757071, 21.0658142, 0.0), 29448268, "service")
     Main.OSMXGraph.Edge(6, 10565, 89783, 1495256391, 8213623939, LLA(52.2756593, 21.0657277, 0.0), LLA(52.2757071, 21.0658142, 0.0), 29448268, "service")



### 3. Building data structures
The data frame for graph metadata is created based on the edges generated in the previous step. Each row
represents a one-way road and includes both node and way IDs from OSM file. Additionally, each node
and way is assigned new IDs (consecutive natural numbers) required for creating graph matrix and spatial
indexing structures. The road graph is represented as a sparse adjacency matrix for efficiency. The sparse
matrix maps ’from id’ and ’to id’ to the road’s id. The ’create road graph’ function also generates KDTree
spatial index for location of the road nodes.


```julia
df = OSMXGraph.edges_to_df(edges)
sparse_index = OSMXGraph.create_sparse_index(df.from_id,df.to_id,df.id)
```


    176262×176262 SparseArrays.SparseMatrixCSC{Int64, Int64} with 358494 stored entries:
    ⎡⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣿⎤
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣟⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢽⣯⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣽⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣾⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣿⣿⣿⣿⡟⣿⣿⣿⣯⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣿⣿⡿⣻⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⡯⣿⣿⣿⣿⡏⢽⣿⣿⣿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣽⣿⣾⣿⣯⣯⡿⣯⣿⣿⣿⣍⣿⣟⣽⣟⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣯⣾⣿⣿⣿⎥
    ⎢⣿⡿⣿⣿⣿⣿⣿⢿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣿⣿⡿⢿⣿⣿⠿⣿⡿⡿⠿⡟⢿⡿⣿⣿⣿⣿⡿⢿⡿⎥
    ⎢⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣾⣿⣿⣿⣷⣿⣿⣿⣾⣿⣿⣿⣿⣿⣷⣷⣿⢿⣾⣿⣿⡿⣿⣾⣿⣿⎥
    ⎣⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿⣿⣿⣿⣿⣷⢿⣿⣿⣿⡷⣿⣿⣿⣿⎦


### 4. Finding Nearest Points
The module uses OSMToolset POI data frame as a base for further calculations. It takes
the POI metadata structure as input and extends it by a new column with the ID of the
nearest road node, using the KDTree spatial index created in the previous step.  


```julia
lats = [node[1].lat for node in values(nodes)]
lons = [node[1].lon for node in values(nodes)]
vals = [node[2] for node in values(nodes)]
road_mtrx = Matrix(transpose([lats lons]))
road_index = OSMXGraph.create_road_index(road_mtrx)
```


    NearestNeighbors.KDTree{StaticArraysCore.SVector{2, Float64}, Distances.Euclidean, Float64, StaticArraysCore.SVector{2, Float64}}
      Number of points: 176262
      Dimensions: 2
      Metric: Distances.Euclidean(0.0)
      Reordered: false


The building process can be shortened if some of the structures are saved as files. The 
module allows for saving and reading both graph metadata and POI data frames as .csv 
files. Additionally, there is an option to save a hashmap with the key as 'OSM node id' 
and the value as [geodetic coordinates, graph ID]. The sparse matrix and KDTree are not 
saved but can be quickly regenerated based on saved files.

## Workflow

The above pipeline is implemented in the functions create_road_graph and add_nearest_road_point. A use case is presented below.


```julia
dir_in = "../data"
road_types = ["motorway", "trunk", "primary", "secondary", 
            "tertiary", "residential", "service", "living_street", 
            "motorway_link", "trunk_link", "primary_link", "secondary_link", 
            "tertiary_link"]  
df, sparse_index, road_index, node_ids = OSMXGraph.create_road_graph("Warszawa.osm", road_types,"Warszawa_graph.csv","Warszawa_nodes.json",dir_in=dir_in)
POI_df = OSMToolset.find_poi(string(dir_in,"/","Warszawa.osm"))
POI_xs = POI_df.lat
POI_ys = POI_df.lon
poi_with_nearest_points = OSMXGraph.add_nearest_road_point(POI_df,POI_xs, POI_ys, road_index, node_ids)
OSMXGraph.save_file(poi_with_nearest_points,"poi_with_nearest_points.csv")
```


    "./poi_with_nearest_points.csv"


## Aknowledgments
This research was funded by National Science Centre, Poland, grant number 2021/41/B/HS4/03349.

<small>The module builds upon two different packages for manipulating OpenStreetMap data: OSMToolset.jl (https://github.com/pszufe/OSMToolset.jl) and OpenStreetMapX.jl (https://github.com/pszufe/OpenStreetMapX.jl). <small>
