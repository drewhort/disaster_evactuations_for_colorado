import geopandas
import matplotlib.pyplot as plt
import momepy
import numpy as np
import networkx as nx
#from contextily import add_basemap
import osmnx as ox
from scipy.interpolate import UnivariateSpline

#potential start points
#71 Blue Mountain Road, Lyons, CO
#start=(40.249325, -105.290504)

#10345 Ute Hwy, Longmont, CO 80504
#stop=(40.20516309783108, -105.10971363190326)

#the code below road calculates capacities based on safe breaking distances and a paper from Moore et al (2013) and then appends this attributed to the graph G.[ code from mike schimidt ]

MOORE_AFTER_BREAK_SPLINE = UnivariateSpline(
    [20, 30, 40, 50, 60, 70, 80, 90, 100],
    [3.9, 6, 11, 18, 27, 39, 54, 58, 84],
)
MOORE_BEFORE_BREAK_SPLINE = UnivariateSpline(
    [20, 30, 40, 50, 60, 70, 80, 90, 100],
    [6, 8, 11, 14, 17, 19, 22, 25, 28],
)

MOORE_SAFE_BREAKING_DISTANCE = lambda x: MOORE_AFTER_BREAK_SPLINE(
    x
) + MOORE_BEFORE_BREAK_SPLINE(x)


def moore(lanes: float, max_speed: float):
    return 1000 * max_speed / MOORE_SAFE_BREAKING_DISTANCE(max_speed) * lanes


#the below function adds capacities to our multidigraph G using the moore method from above

def add_capacities(G, method=moore):
    G = G.copy()
    cap = []
    for u, v, i in G.edges:
        edge_data = G.get_edge_data(u, v, i)
        raw_lanes = edge_data.get("lanes")
        if raw_lanes is None:
            lanes = 1
        elif isinstance(raw_lanes, str):
            lanes = int(raw_lanes) / 2  
        elif isinstance(raw_lanes, list):
            lanes = sum(int(x) for x in raw_lanes) / 2
        edge_data["capacity"] = method(lanes, edge_data["speed_kph"])
        #print('capacity',edge_data["capacity"])
        val = edge_data["capacity"]
        print('val',val)
        #print(edge_data["capacity"])
        cap.append(val)
    return (cap,G)
    

    
#import data from osmnx, can input any city, state, etc.
G = ox.project_graph(ox.graph_from_place('Lyons, Colorado', network_type='drive'))
G_original = ox.project_graph(ox.graph_from_place('Lyons, Colorado', network_type='drive'))
#get rid of intersections that are not actually intersections
G = ox.consolidate_intersections(G, tolerance=10, rebuild_graph=True, dead_ends=True)
G_original = ox.consolidate_intersections(G_original, tolerance=10, rebuild_graph=True, dead_ends=True)
#add edge speeds
G = ox.add_edge_speeds(G)
G_original = ox.add_edge_speeds(G_original)
#add travel times
G = ox.add_edge_travel_times(G)
G_original = ox.add_edge_travel_times(G_original)


#add capacities to G and define edge_data list
edge_data,G_new =add_capacities(G)

#print G as a df
#G_as_df = ox.graph_to_gdfs(G)
#print('gdf',G_as_df)
#print('edge data',edge_data)

#define our weights for the shortest path

#define a start node and a sink node (currently random)
orig, dest = list(G)[10], list(G)[20]

#find initial shortest path using OSMNX from source to sink
shortest_path=ox.shortest_path(G, orig, dest, weight="travel_time")
print('shortest path from s to t', type(shortest_path),shortest_path)

#how to turn a multidigraph into a digraph
#diGraph = ox.utils_graph.get_digraph(G, weight="travel_time")

#shorest path from all nodes to the sink
#shortest_path_all=nx.shortest_path(G, source=None ,target=dest, weight="travel_time")
#print('shortest path from all nodes to t', type(shortest_path_all),shortest_path_all)


#make capacities into an array to append to our edge matrix
capacity=np.floor(np.array(edge_data))
nx.set_edge_attributes(G,capacity,'capacity')
print("edge capacitiy is ", capacity)

#can try to get max flow on a digraph
#flow_value, flow_dict = nx.maximum_flow(diGraph, orig, dest, capacity='capacity')


#make G into an edge matrix
edgematrix=[list(i) for i in G.edges]
edgematrix=np.array(edgematrix)
edgematrix[:,2]=capacity.T[:]
print('edge matrix',type(edgematrix),edgematrix)



def upres(edgematrix, path):

    #calculates flow
    flo = np.amax(edgematrix[:,2])
   
    for i in range(len(path)-1):
        j = path[i]
        k = path[i+1]
        cap = edgematrix[(np.where((edgematrix[:,0]==j) & (edgematrix[:,1]==k))),2]
        print('cap',cap)
        if cap < flo:
            flo = cap
    print('flow',flo)
    
    
    #updates residual graph to give us max flow
    for i in range(len(path)-1):
        j = path[i]
        k = path[i+1]
        #print('edge matrix',np.where((edgematrix[:,0]==j) & (edgematrix[:,1]==k)))
        edgematrix[(np.where((edgematrix[:,0]==j) & (edgematrix[:,1]==k))),2]-=flo
        edgematrix[(np.where((edgematrix[:,0]==k) & (edgematrix[:,1]==j))),2]+=flo
        #print(flo)
        #print('what',i,edgematrix)
    return (edgematrix,flo)

#Initializing the stuff
edgematrix,flo=upres(edgematrix,shortest_path)
print('edge matrix', type(edgematrix), edgematrix)
print('flow',type(flo),flo)
#flo=0

while flo > 0:

    #make the edgematrix into a list
    edgematrix=edgematrix.tolist()
    print('!!!!!edgematrix list!!!!!',type(edgematrix),edgematrix)
    #create an edge list
    edge_list=[l[:2] for l in edgematrix]
    print('!!!edgelist!!!!',type(edge_list),edge_list)
    #empty multidigraph
    G = nx.MultiDiGraph()
    #add edges back in
    G.add_edges_from(edge_list)
    
    #turn edgematrix into an array
    edge_weights=np.array(edgematrix)
    print('!!!!1edge_weights!!!!!', edge_weights)
    #get the weights from the array
    edge_weights = edge_weights[:,2] #this is giving the same capacities from earlier
    print('!!!!2edge_weights!!!!', edge_weights)
    #append weights and capacities to the graph G
    nx.set_edge_attributes(G, edge_weights,'weights')
    #nx.set_edge_attributes(G,capacity,'capacity')
    
    #print('8!!!!!!edge_weights!',edge_weights)

    #recompute the shortest path from s to t in the residual network
    shortest_path=ox.shortest_path(G, orig, dest, weight="travel_time")
    print('9!!!!!!!shorest_path!!',shortest_path)
    fig, ax = ox.plot_graph_route(G_original, shortest_path)
    plt.show()
    
    print('edge matrix!!!!!!!',type(edgematrix),edgematrix)
    #find the max flow of the graph
    edgematrix=np.array(edgematrix)
    edgematrix,flo=upres(edgematrix,shortest_path) #is shortest path really what we want here? 
    print('10!!!!!!!!edge_matrix!',edgematrix)
    #print(nx.graph_to_gdfs(G))
    #G=ox.project_graph(G)
    
    #turn G into an edge matrix
    #edgematrix=[list(i) for i in G.edges]
    #print('1!!!!edgematrix',type(edgematrix),edgematrix)
    #print(edgematrix) #this prints edge matrix, but every weight is 0
    #turn edge matrix into an array
    #edgematrix=np.array(edgematrix)
    #print('2!!!!edgematrix',type(edgematrix),edgematrix)
    #append capacity to edge matrix
    #edgematrix[:,2]=capacity.T[:] #here we need to add new capacities, not the old this is currently  printing old capacities
    #print('3!!!!!',edgematrix)
    
    
    
   
