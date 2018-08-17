# traj-converters
Utilities for converting output of trajectory inference algorithms.

Common formats are described [here](https://github.com/Stuartlab-UCSC/traj-formats)

To use, start off by installing the trajectory inference algorithm you wish to convert. Then clone this repository to your machine at 'some/directory'. Converters use a subset of their trajectory inference algorithm's dependencies, so if an algorithm is installed then the converter should work with no extra hassle. If that's not the case please submit an issue.

## python converters
* scimitar
* wishbone

Set your $PYTHONPATH to point at 'some/directory/traj-converters/src/python'.

Open a python terminal, import the converter and use python's native help(module) for documentation. 

## R converters
* monocle
* slicer
* dpt

The general pattern is to `source("some/directory/traj-converters/src/R/method_convert.r")`

Change *method* to one listed above.

Each *method*_convert.r file has two functions for writing method outputs to a common format.

`write_cell_x_branch(TI_obj, file)` writes a tab delimited matrix to *file* with empty spaces representing NA.

`write_common_json(TI_obj, file)` writes json of the common format to *file*.  

In most cases *TI_obj* is the R object returned by the trajectory inference method. 

Methods that do not return an object require multiple parameters. The parameter names mimic the names of the TI method's function calls used to produce those values. e.g. The TI method slicer has a cell_order() function used for pseudotime assignment. Slicer's converter functions have a parameter *cell_ordering* needed to make the common formats.  

Each *method*_convert.r file also has two functions for creating the R analog data structure.

`to_cell_x_branch(TI_obj)` returns a cell x branch matrix with values representing pseudotime assignment. 

`to_common_list(TI_obj)` returns a list with the schema of the common json format.

