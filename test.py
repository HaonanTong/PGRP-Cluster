import numpy as np
import pandas as pd 
from Visualization import Visualization
from Operations import Operations

ops = Operations()
ops.f_getSignature()
print ops.f_KMEANS_Sig(nClusters=17, isTSNE=True)
print '-----------------------------------------'
print ops.f_KMEANS_Sig(nClusters=17, isTSNE=False)