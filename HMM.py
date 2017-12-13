###################################################################
# HMM: Training, Clustering and Visualization
###################################################################

import numpy as np
import pandas as pd 
from Visualization import Visualization
from Operations import Operations

######################## Training #############################
# Load Data And Prepare For Training
GoI = ['AT2G37090','AT4G18780','AT5G15630','AT5G17420','AT5G44030'] 
ops = Operations()
ops.f_GoI(GoI)


n_components = 3
ops.f_HMM_Discretization(n_components=n_components)
if ops.f_isGoI_same_Z():
	print 'GoI have same pattern'
	print ops.f_GoI_same_Z()
else:
	print 'GoI do not have same pattern'
	raw_input('Press to See State Sequence for GoI...')
	print ops.f_GoI_State()
	raw_input('Adjust Number Of Hidden State To Make Sure GoI Are Clustered Together')


######################## Clustering #############################
X_Prediction, glist = ops.f_Z_Assoc()
# print X_Prediction
print '# Of Gene In The Cluster Assoc With GoI: ' + str(len(X_Prediction))


# Validation With 28 Validated TFs
TP_set, Pred_set = ops.f_HMM_Prediction()
# Intersection For Phisical TFs And Prediction
print 'TPs: ' + str(TP_set)
print 'Predictions: ' + str(Pred_set)

############### Functional Analysis Of Predictions ###############
ops=Operations()
ops.f_get_TSNE()

# Tendency for GoI
ops.f_get_GoI_label()

# Tendency for TP
ops.f_GoI(TP_set)
ops.f_get_GoI_label()

# Tendency for Pred
ops.f_GoI(Pred_set)
ops.f_get_GoI_label()
# raw_input('Press to Cntnue')



############### Visualization ###############
## Plot time-course signature of cluseter
vis = Visualization(Dir = 'Img/HMM')
# Visulaize Trajectory of genes in the cluster assoc with GoI
vis.f_GoI(glist)
vis.f_plot_trajectory(isLegend=False, isGoI=True,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Trajectories_HMM_Cluster.pdf')

# Visulaize Trajectory of Predicted TFs
vis.f_GoI(Pred_set)
vis.f_plot_trajectory(isLegend=True, isGoI=True,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Trajectories_HMM_Pred.pdf')

# Visulaize Validated TFs
vis.f_GoI(TP_set)
vis.f_plot_trajectory(isLegend=True, isGoI=True,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Trajectories_HMM_TP.pdf')


## Plot tSNE signature of cluseter
ops = Operations()
X = ops.f_Inhibitors_upsidedown()
ops = Operations(X)
ops.f_get_TSNE()
vis.f_set_df_tsne(ops.df_tsne)

vis.f_GoI(glist)
vis.f_plot_TSNE(FileName='TSNE_HMM_Cluster.pdf', isTrans=True, isLabel=False, isGoI=True, isOutput=True)
vis.f_GoI(Pred_set)
vis.f_plot_TSNE(FileName='TSEN_HMM_Pred.pdf', isTrans=True, isLabel=False, isGoI=True, isOutput=True)
vis.f_GoI(TP_set)
vis.f_plot_TSNE(FileName='TSEN_HMM_TP.pdf', isTrans=True, isLabel=False, isGoI=True, isOutput=True)





