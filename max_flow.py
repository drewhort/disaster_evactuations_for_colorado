import geopandas
import pandas as pd
import matplotlib.pyplot as plt
import momepy
import numpy as np
import networkx as nx
#from contextily import add_basemap
import osmnx as ox
import matplotlib.cm as cm
from scipy.interpolate import UnivariateSpline

'''
Calculates capacities based on safe breaking distances and a paper from Moore et al (2013) and then appends this attributed to the graph G. [ code adapted from Mike Schimidt ]
'''

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



#this uses the speed_kph attribute of the graph to calculate the moore capacities 
def moore(lanes: float, max_speed: float):
    return 1000 * max_speed / MOORE_SAFE_BREAKING_DISTANCE(max_speed) * lanes


'''
Adds capacities to multidigraph G using moore method from above [ code adapted from Mike Schimidt ]
'''

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
            
        edge_data["capacity"] = int(method(lanes, edge_data["speed_kph"]))
        
        val = edge_data["capacity"]
        #print('val',val)
 
        cap.append(val)
    return (cap,G)

    
'''
MAXIMUM FLOW SHORTEST PATH ALGORITHM (ONE ITERATION) FOR DATA FRAME
'''

def upres(edgematrix, path):

    caps=edgematrix['capacity']
    flo = caps.max()
    #finds the capacities of each edge in the shortest path
    for i in range(len(path)-1):
        j = path[i]
        k = path[i+1]

        #find the rows of the data frame that match the arc (j,k) in the network and the arc (k,j) in the residual network
        forward_index = edgematrix.loc[(edgematrix['source']==j) & (edgematrix['target']==k)]
        backward_index = edgematrix.loc[(edgematrix['source']==k) & (edgematrix['target']==j)]
        
        #If the arcs don't exist in the data frame (i.e. one way road in the original network, then we need to add in the arc to the residual
        if forward_index.empty == True:
             #appends a row of all 0's to the empty dataframe forward_index
             forward_index=forward_index.append(pd.Series(0, index=edgematrix.columns), ignore_index=True)
             #assigns the correct values of j,k, 'capacity' to the arc
             forward_index.at[0,'capacity']= 0
             forward_index.at[0,'source']= j
             forward_index.at[0,'target']= k
             #appends the row to the edgematrix
             edgematrix=pd.concat([edgematrix,forward_index],ignore_index=True)
             #we can now find the index in edgematrix that corresponds to (j,k) so that we can appropriately update capacity after augmenting flow
             forward_index = edgematrix.loc[(edgematrix['source']==j) & (edgematrix['target']==k)]
             
        #this works the same as above, but for edge (k,j)     
        if backward_index.empty == True:
             backward_index=backward_index.append(pd.Series(0, index=edgematrix.columns), ignore_index=True)
             backward_index.at[0,'capacity']= 0
             backward_index.at[0,'source']= k
             backward_index.at[0,'target']= j
             edgematrix=pd.concat([edgematrix,backward_index],ignore_index=True)
             backward_index = edgematrix.loc[(edgematrix['source']==k) & (edgematrix['target']==j)]
        
        #get index for the forward arc from our path
        for row in forward_index.index:
            forward=row

        #get the edge capacity from the graph
        cap = edgematrix.loc[forward].at['capacity']
        
        #if the capacity is lower than the current flow, set flow to be the capacity of the edge
        if cap < flo:
            flo = cap
            
    print('flow',flo)
    
    #augments flow and updates residual graph
    for i in range(len(path)-1):
        j = path[i]
        k = path[i+1]
        
        #again find the location of our indices
        forward_index = edgematrix.loc[(edgematrix['source']==j) & (edgematrix['target']==k)]
        backward_index = edgematrix.loc[(edgematrix['source']==k) & (edgematrix['target']==j)]
        
        #assign the indices 
        for row in forward_index.index:
            forward=row
        
        for row in backward_index.index:
            backward=row 
        
        #augment flow and update the residual network.
        edgematrix.at[forward,'capacity']= edgematrix.loc[forward].at['capacity']-flo
        edgematrix.at[backward,'capacity']= edgematrix.loc[backward].at['capacity']+flo

        
    return (edgematrix,flo)    
    

'''
Building the Graph
'''
    
#import data from osmnx, can input any city, state, etc.
G = ox.project_graph(ox.graph_from_place('Butte, California', network_type='drive'))

#get rid of intersections that are not actually intersections
G = ox.consolidate_intersections(G, tolerance=10, rebuild_graph=True, dead_ends=True)

#add edge speeds
G = ox.speed.add_edge_speeds(G)

#add travel times
G = ox.speed.add_edge_travel_times(G)

#add capacities (computed using moore method)
edge_data,G =add_capacities(G)

#G changes, so we want to have an original network to later plot the shortest paths on
G_original=G.copy()



'''
Plotting the Graph with edge colors corresponding to edge attributes
'''

#plot G with capacities as colors for edges
#ec2=ox.plot.get_edge_colors_by_attr(G, attr="capacity")
#fig, ax = ox.plot_graph(G,edge_color=ec2, node_size=2)
#plt.show()


#plot G with the edges colored by speed kph
ec = ox.plot.get_edge_colors_by_attr(G_original, attr="speed_kph")
#fig, ax = ox.plot_graph(G_original, edge_color=ec,node_size=2)
#plt.show()


#plot G with the edges colored by travel time
#ec = ox.plot.get_edge_colors_by_attr(G_original, attr="travel_time")
#fig, ax = ox.plot_graph(G_original, edge_color=ec,node_size=2)
#plt.show()
   

'''
Find the inital flow (Shortest Path)
You need to choose a orig and dest, and also make sure that the respective area is put into line 148.
'''

#define a start node and a sink node

#boulder county (use boulder county as map)
#this origin is where the marshall fire started
#orig=ox.distance.nearest_nodes(G,479800,4424260,return_dist=False)
#destination is the evacuation center listed for marshall fire
#dest=ox.distance.nearest_nodes(G,489130,4429760,return_dist=False)
#alternative origin that is in superior which was higly affected by the marshall fire
#orig=ox.distance.nearest_nodes(G,485830,4423270,return_dist=False)

#pueblo west (use pueblo county as the map)
#not the right origin/dest I think
#orig=ox.distance.nearest_nodes(G,5255180,4250280,return_dist=False)
#orig=ox.distance.nearest_nodes(G,514380,4242010,return_dist=False)

#Another Pueblo West
#orig=ox.distance.nearest_nodes(G,514930, 4241970,return_dist=False)
#dest=ox.distance.nearest_nodes(G,530020,4232490,return_dist=False)

#Lyons, CO (run on boulder county) 
#orig=ox.distance.nearest_nodes(G,477230,4453570,return_dist=False)
#dest=ox.distance.nearest_nodes(G,476140,4430890,return_dist=False)

#Ken Caryl (Jefferson County, CO)
#orig=ox.distance.nearest_nodes(G,489685,4379999,return_dist=False)
#dest=ox.distance.nearest_nodes(G,494940,4393420,return_dist=False)

#Heildsburg, CA (run on Sonoma County, CA)
#orig=ox.distance.nearest_nodes(G,512007,4273277,return_dist=False)
#santa rosa destination
#dest=ox.distance.nearest_nodes(G,525300,4255270,return_dist=False)
#petaluma destination
#dest=ox.distance.nearest_nodes(G,532060,4233990,return_dist=False)

#Paradise (Butte County, CA)
#where the camp fire started
orig=ox.distance.nearest_nodes(G,621540,4406520,return_dist=False)
#chico was one of the places you could evacuate to
dest=ox.distance.nearest_nodes(G,599230,4400090,return_dist=False)


#define empty list to store shortest paths
shortest_path_list=[]

#find initial shortest path using OSMNX from source to sink
shortest_path=ox.distance.shortest_path(G, orig, dest, weight="travel_time")
print('shortest path from s to t', type(shortest_path),shortest_path)

#add the initial shortest path to the list
shortest_path_list.append(shortest_path)

#plot the 1st shortest path
#fig, ax = ox.plot_graph_route(G_original, shortest_path,edge_color=ec, node_size=2)
#plt.show()

#transform graph into a dataframe
edgematrix = nx.to_pandas_edgelist(G)

#egdgematrix=ox.utils_graph.graph_to_gdfs(G,)



'''
Maximum Flow Algorithm
'''
#Initializing the original flow using the first shortest path found above
edgematrix,flo=upres(edgematrix,shortest_path)

#i is the counter 
i=1

#while loop iterates through each shortest augmenting path
while shortest_path is not None:
    i+=1
    
    #convert dataframe back to graph to get the shortest path using the osmnx shortest path function
    G = nx.from_pandas_edgelist(edgematrix,'source','target',edge_attr='travel_time') #changed from cap to travel_time not sure if this makes a difference
    G_shorty=G
    
    #want to temporarily delete 0 capacity edges before we send to the shortest path
    zeros_index = edgematrix.loc[(edgematrix['capacity']==0.0)]
    for j in zeros_index.index:
        #print(G_shorty[edgematrix.iloc[j,0]][edgematrix.iloc[j,1]]['capacity'])
        G_shorty.remove_edge(edgematrix.iloc[j,0],edgematrix.iloc[j,1])
   

    #recompute the shortest path from s to t in the residual network with weights travel_time
    shortest_path=ox.distance.shortest_path(G_shorty, orig, dest, weight='travel_time')
    
    #Append non-empty shortest paths to the shortest path list
    if shortest_path is not None:
        shortest_path_list.append(shortest_path)
        print('iteration',i,'shortest path', shortest_path)
        
        #plot each shortest path
        #fig, ax = ox.plot_graph_route(G, shortest_path,edge_color=ec, node_size=2)
        #plt.show()
       
        edgematrix,flo=upres(edgematrix,shortest_path) 

    
#print all of our different paths.   
print('shortest path list',shortest_path_list, len(shortest_path_list))

'''
We want to only plot the paths that exist in the original network
'''
sp_list=[]
for path in shortest_path_list:
    sp=[]
    #print('path',path)
    for i in range(len(path)-1):
        j=path[i]
        k=path[i+1]
        sp.append((j,k))
    #print('sp',sp)    
    t=all(G_original.has_edge(*l) for l in sp)
    #print(t)
    if t==True:
        sp_list.append(path)
    
print(sp_list)

'''
Get the maximum flow value
'''
max_flow=edgematrix.loc[orig].at['capacity']

print('max_flow',max_flow)


'''
Gets n evenly spaced colors to map the different path with from the colormap 'hsv'
'''

def get_colors(n, cmap='hsv', start=0., stop=1., alpha=1.):

    colors = [cm.get_cmap(cmap)(x) for x in np.linspace(start, stop, n)]
    colors = [(r, g, b, alpha) for r, g, b, _ in colors]
    return colors


'''
Plot all the shortest paths
'''
if len(sp_list)>1:
    cl=get_colors(len(sp_list))
    fig, ax = ox.plot_graph_routes(G_original, sp_list, route_colors=cl, route_linewidth=6, node_size=2,edge_color=ec)
    plt.show()
else:  
    fig, ax = ox.plot_graph_route(G_original, sp_list[0], route_color='y', route_linewidth=6, node_size=2,edge_color=ec)
    plt.show()
    

    
