###################################################################
# Part1: KMeans on Time-course signature
###################################################################

import numpy as np
import pandas as pd 
from Visualization import Visualization
from Operations import Operations
## Result
ops = Operations()
ops.f_getSignature()
ops.f_KMEANS_Sig(nClusters=17, isTSNE=False)

List_Cluster, List_Pred, List_TP = ops.f_KMeans_Result()
print '# of Genes in the cluster Assoc With GoI: ' + str(len(List_Cluster))
print List_Cluster
print '# of Predictions: ' + str(len(List_Pred))
print List_Pred
print '# of True Positive: ' + str(len(List_TP))
print List_TP

## Visualization by trajectories
vis = Visualization(Dir = 'Img/tsig')
vis.f_GoI(List_Cluster)
vis.f_plot_trajectory(isLegend=False, isGoI=True,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Trajectories_KMeans_tsig_cluster.pdf')

vis.f_GoI(List_Pred)
vis.f_plot_trajectory(isLegend=True, isGoI=True,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Trajectories_KMeans_tsig_Pred.pdf')

vis.f_GoI(List_TP)
vis.f_plot_trajectory(isLegend=True, isGoI=True,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Trajectories_KMeans_tsig_TP.pdf')

## Visualization by tSNE-signature
vis.f_set_df_tsne(ops.sig_tSNE)

vis.f_GoI(List_Cluster)
vis.f_plot_TSNE(FileName='TSNE_KMeans_tsig_Cluster.pdf', isTrans=True, isLabel=False, isGoI=True, isOutput=True)
vis.f_GoI(List_Pred)
vis.f_plot_TSNE(FileName='TSEN_KMeans_tsig_Pred.pdf', isTrans=True, isLabel=False, isGoI=True, isOutput=True)
vis.f_GoI(List_TP)
vis.f_plot_TSNE(FileName='TSEN_KMeans_tsig_TP.pdf', isTrans=True, isLabel=False, isGoI=True, isOutput=True)


vis.f_setlabel(ops.Result_Kmeans)
vis.f_plot_TSNE(FileName='TSEN_KMeans_tsig_17C.pdf', isTrans=False, isLabel=True, isGoI=False, isOutput=True)

## Functional Analysis
ops.f_get_TSNE()
# Tendency for GoI
print 'GoI label: \n' + ops.f_get_GoI_label().to_string()
# Tendency for TP
# ops.f_GoI(List_TP)
# print ops.f_get_GoI_label().to_string()
# Tendency for Pred
ops.f_GoI(List_Pred)
print 'Predictions label: \n' + ops.f_get_GoI_label().to_string()



