## Currently there is only one converter from the output of the scimtar pipeline to the v1 of the common trajectory format.

Here's some example usage. If your interested in using this but don't know how to get "scimatar pipeline" output, please open an issue.  
```
#The module must be in your python path.
module_dir='path/to/src/python'
#import sys
#sys.path.append(module_dir)

import pickle
import scimitar2json
import os

data_path = os.path.join(
    module_dir,
    '../../examples/data/RG_PC_pickle_test.pi'
)

# Grab the data in the example dir.
scimitar_model, diff_map = pickle.load(open(data_path, 'rb'))

# Make a python dictionary of the common format.
common_dict = scimitar2json.make_dict(scimitar_model, diff_map)
# Save the dictionary to json format.
scimitar2json.save_to_json(common_dict, data_path)
```
