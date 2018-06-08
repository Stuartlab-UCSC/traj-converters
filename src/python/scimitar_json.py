"""
Contact info: lseninge@ucsc.edu 

Script converting SCIMITAR output from SC-SUITE pipeline to consensus
trajectory format.

See help(scimitar_jsona) for extensive documentation.

Example usage:

#The module must be in your python path.
module_dir='path/to/src/python'
#import sys
#sys.path.append(module_dir)

import pickle
import scimitar_json
import os

data_path=os.path.join(module_dir, '../../examples/data/RG_PC_pickle_test')

scimitar_model, diff_map=pickle.load(open(data_path, 'rb'))

# Make a python dictionary of the common format.
common_dict = scimitar_json.make_dict(scimitar_model, diff_map)
# Save the dictionary to json format.
scimitar_json.save_to_json(common_dict,data_path)
"""

import numpy as np
import networkx as nx
import json
from numpy.linalg import norm


def make_dict(scimitar_output=None,diff_map=None):
    """
    This fuction convert SCIMITAR output and probability assignment to consensus JSON format
    :param scimitar_output: dictionary containing SCIMITAR model, probabilities and cell assignment to node
    :param diff_map: diffusion map coordinates of cells

    :return common_dico: A dictionary containing consensus fields for trajectory description
    """
    common_dico={'nodes':{'nodeId':_nodeId(scimitar_output),'dimReduct':_dimReduct(scimitar_output).tolist()},
                 'edges':{'edgeId':_edgeId(scimitar_output)[0],'nodeId1':_edgeId(scimitar_output)[1],'nodeId2':_edgeId(scimitar_output)[2],'direction':_direction(scimitar_output)},
                 'cellMapping':{'cellId':_cellMapping(scimitar_output,diff_map)[0],'branchId':_cellMapping(scimitar_output,diff_map)[1],'psuedotime':_psuedotime(scimitar_output,diff_map),'curveDistance':_cellMapping(scimitar_output,diff_map)[2]}
                }
    return common_dico


def save_to_json(dico, path):
    """
    This function save the dictionary output by the make_dict function to JSON format
    :param dico: dictionary output by the make_dict function
    :param path: /path/to/save/file
    
    :return NA
    """
    #Encoder
    class MyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            else:
                return super(MyEncoder, self).default(obj)

    #Save
    json_format = json.dumps(dico, cls=MyEncoder, indent=4)
    with open(path, 'w') as outfile:
        json.dump(json_format, outfile)


def _nodeId(scimitar_output):
    """
    This function get node ID from the SCIMITAR output inside the SC-SUITE pipeline
    :param scimitar_output: dictionary containing SCIMITAR model, probabilities and cell assignment to node
    
    :return A list of node ID corresponding to cell assignment to node
    """
    nodeId=list(set(np.unique(scimitar_output[2]).tolist()))
    return nodeId

def _dimReduct(scimitar_output):
    """
    This function gets node positions in the reduced space (Diffusion Map)
    :param scimitar_output: dictionary containing SCIMITAR model, probabilities and cell assignment to node
    
    :return A list of 2-element lists containing x,y coordinates of each node in the reduced space 
    """
    dimReduct=scimitar_output[0].__dict__['node_positions']
    return dimReduct

def _edgeId(scimitar_output):
    """
    This function gets edges of the branching tree output by SCIMITAR and its associated nodes
    :param scimitar_output: dictionary containing SCIMITAR model, probabilities and cell assignment to node
    
    :return edgeId:A list of tuple corresponding to edges (e.g: (1,2))
    :return NodeId1: list of first node of edge
    :return NodeId2: list of second node of edge
    """
    edgeId=[e for e in scimitar_output[0].__dict__['graph'].edges]
    NodeId1=[k[0] for k in edgeId]
    NodeId2=[k[1] for k in edgeId]
    return edgeId, NodeId1, NodeId2

def _direction(scimitar_output):
    """
    There is no direction to the edges in SCIMITAR model yet
    """
    direction=['undirected']*len(_edgeId(scimitar_output)[0])
    return direction

def _cellMapping(scimitar_output,diff_map):
    """ This function gets each cell name and assign it to it's closest edge in the branching tree.
    SCIMITAR assign each cell to a node by default. To assign to closest edge, we use node coordinates and dot product
    to get orthogonal distances to edges and assign each cell to its closest edge
    :param scimitar_output: dictionary containing SCIMITAR model, probabilities and cell assignment to node
    :param diff_map: diffusion map coordinates of cells
    
    :return cellId: cell name in the original matrix
    """
    #cellId
    cellId=diff_map.index.tolist()

    #branchId
    #Get distance of cells to trajectory branches
    array_dist=np.zeros(shape=(len(diff_map),len(_edgeId(scimitar_output)[0])))
    i=0
    #Iterate on each cell
    for index, row in diff_map.iterrows():
        p3=[row[0],row[1]]
        dist_list=[]
        #Find distance of each cell to each edge (orthogonal projection)
        for edge in _edgeId(scimitar_output)[0]:
            node1,node2=edge
            p1, p2 = _dimReduct(scimitar_output)[node1], _dimReduct(scimitar_output)[node2][1]
            d = norm(np.cross(p2-p1, p1-p3))/norm(p2-p1)
            dist_list.append(d)
        array_dist[i]=dist_list
        i+=1
    #Assign cell to its closest branch + save distance to this branch
    branchId=[]
    curveDistance=[]
    for i in range(array_dist.shape[0]):    
            branchId.append(_edgeId(scimitar_output)[0][int(np.argmin(array_dist[i, :]))])
            curveDistance.append(np.min(array_dist[i, :]).item())
            
    return cellId, branchId, curveDistance

def _psuedotime(scimitar_output, diff_map):
    """
    There is no pseudotime assignment to cells in Branching scimitar yet.
    """
    psuedotime=['NA']*len(diff_map)
    return psuedotime



