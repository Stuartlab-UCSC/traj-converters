"""
Module that converts wishbone object to common formats.

As of this writing two common formats have been established: json,
    and a cell x branch matrix with values representing pseudotime
    as a tab delimited matrix.

Descriptions if the common formats can be found here:
    https://github.com/Stuartlab-UCSC/traj-formats

Star topology is assumed for the graph to minimize heuristics.

The output_* functions in this file will write the common format to a
    designated file from a wishbone object.

The to_* functions will convert the wishbone object to a python data
    structure (dict, pandas dataframe)

Notebook tutorial for constructing a wishbone object can be found here:
    http://nbviewer.jupyter.org/github/ManuSetty/wishbone/blob/master/notebooks/Wishbone_for_single_cell_RNAseq.ipynb

Our examples start from the In [21] of that notebook.
>>> wb = wishbone.wb.Wishbone(scdata)
>>> wb.run_wishbone(start_cell='W30258', components_list=[1, 2], num_waypoints=150)

# Output the cell x branch matrix to a file:
>>> output_cell_x_branch(wb, "path/to/output.tab")

# Output the common json format to a file:
>>> output_common_json(wb, "path/to/output.json")

# Make the cell x branch pandas dataframe:
>>> cellXbranch = to_cell_x_branch(wb)

# Make a dictionary representing the common format:
>>> commonDict = to_common_dict(wb)
"""
"""

"""
import pandas as pd
import networkx as nx
import numpy as np
import json


def write_common_json(wishbone_obj, filepath):
    """
    Writes common json format to file.
    :param wishbone_obj: A wishbone object.
    :param file: Path to json output.
    :return: None
    """
    common_dict = to_common_dict(wishbone_obj)
    with open(filepath, 'w') as fp:
        json.dump(common_dict, fp)


def write_cell_x_branch(wishbone_obj, filepath):
    """
    Writes the cell_x_branch format.
    :param wishbone_obj: A wishbone object.
    :param filename: Path tab delimited matrix output.
    :return: None
    """
    return to_cell_x_branch(wishbone_obj).to_csv(
        filepath,
        sep="\t"
    )


def to_cell_x_branch(wishbone_obj):
    """
    Make a cell x branch matrix from a wishbone object.
    :param wishbone_obj: A wishbone object.
    :return: dataframe cell x branch with values of pseudotime.
    """
    branches = wishbone_obj._branch
    pseudotime = wishbone_obj._trajectory
    series = []
    for branch in set(branches):
        cells_on_branch = branches.index[branches == branch]
        branch_name = "branch_" + str(branch)

        series.append(
            pd.Series(pseudotime[cells_on_branch], name=branch_name)
        )

    return pd.concat(series, axis=1)


def to_common_dict(wishbone_obj):
    """
    Create python dictionary reprsenting the common json format.
    :param wishbone_obj: A wishbone object.
    :return: (dict) common format dictionary.
    """
    branches = wishbone_obj._branch
    pseudotime = wishbone_obj._trajectory

    common_dict = {}
    graph = make_graph(wishbone_obj)
    common_dict["nodes"] = {"nodeId": list(graph.nodes())}

    edge_info = [(e[0], e[1], edge_id(e[0], e[1])) for e in graph.edges()]
    nodeId1, nodeId2, edgeId = zip(*edge_info)
    common_dict["edges"] = {
        "edgeId": edgeId,
        "nodeId1": nodeId1,
        "nodeId2": nodeId2
    }

    cellIdCM = []
    edgeIdCM = []
    pseudotimeCM = []
    for edge in graph.edges():
        node1, node2 = edge
        edgeId = edge_id(node1, node2)
        #TODO: less fragile, this implementation cant stand if the
        # 'stem' node has a '_' in it, and the end node must have
        # an '_' in it.
        branch_number = get_branch_number([e for e in edge if "_" in e][1])
        cells_on_branch = branches.index[branches == branch_number]
        cellIdCM.extend(cells_on_branch)
        pseudotimeCM.extend(pseudotime[cells_on_branch])
        edgeIdCM.extend(
            np.repeat(edgeId, len(cells_on_branch))
        )

    common_dict["cellNapping"] = {
        "edgeId": edgeIdCM,
        "cellId": cellIdCM,
        "pseudotime": pseudotimeCM
    }

    return common_dict


def make_graph(wishbone_obj):
    """
    Makes the simplified network x object representing the topology.
    :param wishbone_obj: A wishbone object.
    :return: (networkx Graph)
    """
    n_branches = len(set(wishbone_obj._branch))
    graph = nx.star_graph(n_branches)
    name_map = dict(
        zip(
            range(4),
            # if you remove the '_' it'll break to_common_dict
            map(lambda x: "end_" + str(x), graph.nodes())
        )
    )
    name_map[0] = "stem"
    graph = nx.relabel_nodes(graph, name_map)
    return graph


def edge_id(node1, node2):
    """
    Creates an edge id from two node ids.
    :param node1: (string) node id.
    :param node2: (string) node id.
    :return: (string) edge id.
    """
    return node1 + "_" + node2


def get_branch_number(nodeId):
    """
    Parses a node Id and graphs the branch number.
    TODO: Will not work if you pass it the 'stem node'.
    :param nodeId: (string) node id.
    :return: (int) branch number.
    """
    return int(nodeId.split("_")[-1])