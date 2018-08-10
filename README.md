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

Each *method*_convert.r file has two methods for writing method outputs to a common format.

`write_cell_x_branch(...)`

`write_common_json(...)`
