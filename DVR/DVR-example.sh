#! /bin/bash

# Once the proper environment variables are set, and once DVR is running,
# this script provides a basic example of how to perform actions in DVR

dvrcomm 'load "/path/to/data/2d_example.sdf" > example'
dvrcomm 'example_c = coarsen(example)'
dvrcomm 'save example_c > "/tmp/example_c.sdf"'
dvrcomm 'exit'
