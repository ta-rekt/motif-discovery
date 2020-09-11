import numpy as np
import time
import pickle

import functions_regular as fnr
import functions_subtrees as fns

# this test script is going to compare the runtime of fnr and fns on a
# test network from the community_fitNet corpus. it will write the output
# to a file named test.out. the script will run on a single core.
