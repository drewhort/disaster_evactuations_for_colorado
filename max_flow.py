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


I think something has gone wrong here with there aren't indices for when we are going backward on edges, because they don't exist in the DF. Need to add them in somehow.
'''

def upres(edgematrix, path):
    #print('!!!!!edgematrix',edgematrix)
    #edgematrix=edgematrix.copy()
    caps=edgematrix['capacity']
    #print('caps',caps)
    flo = caps.max()
    #print('capacity flow', flo)
    #print(edgematrix.columns)
    #print('!!!!!path',path)
    for i in range(len(path)-1):
        #print(path)
        j = path[i]
        k = path[i+1]

        
        forward_index = edgematrix.loc[(edgematrix['source']==j) & (edgematrix['target']==k)]
        backward_index = edgematrix.loc[(edgematrix['source']==k) & (edgematrix['target']==j)]
        
        if forward_index.empty == True:
             #print('before before!!',forward_index)
             forward_index=forward_index.append(pd.Series(0, index=edgematrix.columns), ignore_index=True)
             forward_index.at[0,'capacity']= 0
             forward_index.at[0,'source']= j
             forward_index.at[0,'target']= k
             #print('before forward!!!',forward_index)
             edgematrix=pd.concat([edgematrix,forward_index],ignore_index=True)
             forward_index = edgematrix.loc[(edgematrix['source']==j) & (edgematrix['target']==k)]
             #print('after forward!!!',forward_index)
             
             
        if backward_index.empty == True:
             backward_index=backward_index.append(pd.Series(0, index=edgematrix.columns), ignore_index=True)
             #print('before before!! backward',backward_index)
             backward_index.at[0,'capacity']= 0
             backward_index.at[0,'source']= k
             backward_index.at[0,'target']= j
             #print('before!!! backward',backward_index)
             edgematrix=pd.concat([edgematrix,backward_index],ignore_index=True)
             backward_index = edgematrix.loc[(edgematrix['source']==k) & (edgematrix['target']==j)]
             #print('after!!backward',backward_index)
                    
        #print('!!!!!! forward index', forward_index)
        for row in forward_index.index:
            forward=row
            #print('forward',forward)
        #print('!!!!!! backward index', backward_index)
        for row in backward_index.index:
            backward=row
            #print('backward',backward)
        cap = edgematrix.loc[forward].at['capacity']
       
        #print('cap',type(cap),cap)
        
        if cap < flo:
            flo = cap
            
    print('actual flow',type(flo),flo)
    
    #updates residual graph to give us max flow
    for i in range(len(path)-1):
        j = path[i]
        k = path[i+1]
        
        #forward=0
        #backward=0
        
        forward_index = edgematrix.loc[(edgematrix['source']==j) & (edgematrix['target']==k)]
        backward_index = edgematrix.loc[(edgematrix['source']==k) & (edgematrix['target']==j)]
        
        for row in forward_index.index:
            forward=row
        
        for row in backward_index.index:
            backward=row
        
        #print(edgematrix.loc[forward].at['capacity'])
        edgematrix.at[forward,'capacity']= edgematrix.loc[forward].at['capacity']-flo
        edgematrix.at[backward,'capacity']= edgematrix.loc[backward].at['capacity']+flo
        #print(edgematrix.loc[forward].at['capacity'])
        
    return (edgematrix,flo)    
    

'''
Building the Graph
'''
    
#import data from osmnx, can input any city, state, etc.
G = ox.project_graph(ox.graph_from_place('Santa Rosa, California', network_type='drive'))

#get rid of intersections that are not actually intersections
G = ox.consolidate_intersections(G, tolerance=10, rebuild_graph=True, dead_ends=True)

#add edge speeds
G = ox.speed.add_edge_speeds(G)

#add travel times
G = ox.speed.add_edge_travel_times(G)

#add capacities (computed using moore method)
edge_data,G =add_capacities(G)
G_original=G.copy()

G_orig_nodes_df,G_orig_edges_df=ox.utils_graph.graph_to_gdfs(G)


'''
Plotting the Graph with edge colors corresponding to edge attributes
'''

#plot G with capacities as colors for edges
#ec2=ox.plot.get_edge_colors_by_attr(G, attr="capacity")
#fig, ax = ox.plot_graph(G,edge_color=ec2, node_size=2)
#plt.show()


#plot G with the edges colored by travel time
#ec = ox.plot.get_edge_colors_by_attr(G_original, attr="speed_kph")
#fig, ax = ox.plot_graph(G_original, edge_color=ec,node_size=2)
#plt.show()


#plot G with the edges colored by travel time
ec = ox.plot.get_edge_colors_by_attr(G_original, attr="travel_time")
#fig, ax = ox.plot_graph(G_original, edge_color=ec,node_size=2)
#plt.show()
   

'''
Find the inital flow (Shortest Path)
'''

#define a start node and a sink node (currently random)
#Lyons random ex
orig, dest = list(G)[8], list(G)[15]
#Boulder County random ex
#orig, dest = list(G)[1], list(G)[9533]

#paradise only
#orig,dest=list(G)[530],list(G)[531]
#butte county
#orig,dest=list(G)[5031],list(G)[995] # also works for boulder county

#boulder county, start is near where marshall fire was and sink is near where and evacuation center was
#orig=ox.distance.nearest_nodes(G,486600,4421340,return_dist=False)
#dest=ox.distance.nearest_nodes(G,482360,4433360,return_dist=False)

#pueblo west (use pueblo county as the map)
#not the right origin/dest I think
#orig=ox.distance.nearest_nodes(G,5255180,4250280,return_dist=False)
#orig=ox.distance.nearest_nodes(G,514380,4242010,return_dist=False)

#good pueblo west ones!!!
#orig=ox.distance.nearest_nodes(G,514930, 4241970,return_dist=False)
#dest=ox.distance.nearest_nodes(G,530020,4232490,return_dist=False)

#Lyons, CO (run on boulder county) !!!! so good!!
#orig=ox.distance.nearest_nodes(G,477230,4453570,return_dist=False)
#dest=ox.distance.nearest_nodes(G,476140,4430890,return_dist=False)


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


#Initializing the original flow using the first shortest path found above
edgematrix,flo=upres(edgematrix,shortest_path)


'''
While loop iterates through each augmentation of the max flow
'''


i=1
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
   
   
    if shortest_path is not None:
        shortest_path_list.append(shortest_path)
        print('iteration',i,'shortest path', shortest_path)
        
        #plot each shortest path
        #fig, ax = ox.plot_graph_route(G, shortest_path,edge_color=ec, node_size=2)
        #plt.show()
       
        edgematrix,flo=upres(edgematrix,shortest_path) 

    
 #print all of our different paths.   
print('shortest path list',shortest_path_list, type(shortest_path_list), len(shortest_path_list))

#need to get all the paths that exist in original network G
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
        
#shortest_path_list.pop()
print(sp_list)

max_flow=edgematrix.loc[orig].at['capacity']

print('max_flow',max_flow)


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
    
#G=ox.utils_graph.graph_from_gdfs(G_orig_nodes_df,edgematrix,graph_attrs=G_original.crs)   
    
#cl=get_colors(len(shortest_path_list))
#fig, ax = ox.plot_graph_routes(G, shortest_path_list, route_colors=cl, route_linewidth=6, node_size=2,edge_color=ec)
#plt.show()    
    
