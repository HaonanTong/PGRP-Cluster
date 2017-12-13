import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from ggplot import *
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import itertools
import pandas as pd
import random
from hmmlearn import hmm


random.seed(0)
np.random.seed(0)


class Operations(object):
	def __init__(self, df=pd.DataFrame()):
		if not df.empty:
			self.df = df # Normalized Expr
		else:
			db = pd.read_csv('myDB_Validation.csv',index_col = 'id')
			# Contains informations of 2421 DEGs
			X1 = db.loc[:,'t0':'t6']
			X2 = db.loc[:,'TF':'Chen_isTar']
			self.GoI = ['AT2G37090','AT4G18780','AT5G15630','AT5G17420','AT5G44030'] 

			X_GoI = X1.loc[self.GoI,'t0':'t6']
			X_pTFs = X1.loc[X2.TF==1,:]
			X = pd.concat([X_GoI, X_pTFs])
			self.df = X

			ValSet = pd.read_csv('NanValidationSet.csv', index_col = 'Transcription factor')
			# Contains differentially expressed causal relationahip for
			# 8 nTFs correspond to cell wall organization:
			#'AT2G22420', 'AT2G28110', 'AT2G37090', 'AT3G18660'
			#'AT4G18780', 'AT5G15630', 'AT5G17420', 'AT5G44030'
			TFs_List = []
			for i in range(0,len(ValSet)):
				if ValSet.Promoter[i] in self.GoI:
					TFs_List.append(ValSet.index[i])
			TFs_Set = set(TFs_List)
			TFs_newList = list(TFs_Set) # Another 28 TFs that is Genes of Interest
			self.TF_Val = TFs_newList

		# self.GoI = ['AT2G37090','AT4G18780','AT5G15630','AT5G17420','AT5G44030'] 
		self.n_gene = len(self.df)
		self.n_time = len(self.df.columns)
		self.df_tsne = []
		self.Result_Kmeans = []
		self.label = []
		self.GoI_cluster_index = []
		self.HMM_Z = pd.DataFrame()
		self.sig_timecourse = pd.DataFrame()
		self.sig_tSNE = pd.DataFrame()

	def f_getSignature(self):
		self.sig_timecourse = self.f_Inhibitors_upsidedown()
		
		tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=5000,random_state=0)
		tsne_results = tsne.fit_transform(self.sig_timecourse)
		self.sig_tSNE = self.df.copy()
		self.sig_tSNE['x-tsne'] = tsne_results[:,0]
		self.sig_tSNE['y-tsne'] = tsne_results[:,1]

	def f_GoI(self,List_GoI):
		self.GoI = List_GoI

	def f_setDF(self,X):
		if len(X)==self.n_gene:
			self.df=X
		else:
			raw_input('Caution! number of gene is not consistent')
			self.n_gene=len(X)

	def f_Inhibitors_upsidedown(self, isSummary=False):
		self.f_get_TSNE()
		self.f_get_KMEANS(nClusters=2, isTSNE=True)
		X = self.f_normalization()

		label = self.Result_Kmeans
		X['label']=label
		X_upper = -1 * X.loc[X.label==0,'t0':'t6']
		X_lower = X.loc[X.label==1,'t0':'t6']
		if isSummary:
			print '# if potential inhibitors: ' + str(len(X_upper))
			print '# if potential acivators: ' + str(len(X_lower) - len(self.GoI))
			raw_input('Press to Continue...')
		X = pd.concat([X_upper, X_lower])
		X = X.loc[self.df.index,:]
		return X

	def f_normalization(self,isTime=True,isObj=False):
		if isTime:
			dfs = self.df.T
			dfs = ( dfs - dfs.mean() ) / dfs.std()
			return dfs.T

		if isObj:
			dfs = self.df
			dfs = ( dfs - dfs.mean() ) / dfs.std()
			return dfs

	def f_get_TSNE(self):
		df = self.f_normalization(isTime = True)
		tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=5000,random_state=0)
		tsne_results = tsne.fit_transform(df)
		self.df_tsne = df.copy()
		self.df_tsne['x-tsne'] = tsne_results[:,0]
		self.df_tsne['y-tsne'] = tsne_results[:,1]
		# print self.df_tsne


	def f_get_GoI_label(self):
		# effects of GoI, activators or repressors
		self.f_get_KMEANS(isTSNE=True)
		X = self.df
		X['label'] = self.Result_Kmeans
		X_GoI = X.loc[self.GoI,:]
		return X_GoI.label

	def f_get_KMEANS(self, nClusters=2, isTSNE=True):
		if isTSNE:
			tmp = self.df_tsne.loc[:,['x-tsne','y-tsne']]
			kmeans = KMeans(n_clusters=nClusters, random_state=0).fit(tmp)
			self.Result_Kmeans = kmeans.labels_

		if not isTSNE:
			kmeans = KMeans(n_clusters=nClusters, random_state=0).fit(self.df)
			self.Result_Kmeans = kmeans.labels_

		return kmeans.labels_

	def f_KMEANS_Sig(self, nClusters=2, isTSNE=True):
		if isTSNE:
			tmp = self.sig_tSNE.loc[:,['x-tsne','y-tsne']]
			kmeans = KMeans(n_clusters=nClusters, random_state=0).fit(tmp)
			self.Result_Kmeans = kmeans.labels_

		if not isTSNE:
			tmp = self.sig_timecourse
			kmeans = KMeans(n_clusters=nClusters, random_state=0).fit(tmp)
			self.Result_Kmeans = kmeans.labels_

		return kmeans.labels_

	def f_grep_GoI(self,USR_DEF_label = [],cluster_GoI = []):
		X = self.df.copy()
		X['label'] = USR_DEF_label
		X_GoI = X.loc[X.label==cluster_GoI,'t0':'t6']
		self.GoI = X_GoI.index
		return self.GoI

	def f_query_where_GoI(self, whichMethod='KMeans'):
		X = (self.df).copy()
		if whichMethod=='KMeans':	
			X['label'] = self.Result_Kmeans
			X_GoI = X.loc[self.GoI,:]
			self.GoI_cluster_index = X_GoI.label.unique()
			return X_GoI

	def f_KMeans_Result(self):
		if not self.GoI_cluster_index:
			self.f_query_where_GoI()

		X = self.df.copy()
		X['label'] = self.Result_Kmeans
		X_Cluster = X.loc[X.label==self.GoI_cluster_index[0],:]
		List_Cluster = X_Cluster.index.tolist()
		List_Pred = list(set(List_Cluster) - set(self.GoI))
		List_TP = list(set(List_Pred)&set(self.TF_Val) )

		return List_Cluster, List_Pred, List_TP

	def f_HMM_Discretization(self, n_components=2):
		# mydf = ops.f_normalization(isTime=True,isObj=False)
		mydf = self.f_Inhibitors_upsidedown()
		# ---- Prepare For Training HMM Model ----
		X = []
		for i in range(0,self.n_gene):
			for j in range(0,self.n_time):
				X.append([mydf.iloc[i,j]])
		lengths = self.n_time * np.ones(self.n_gene)
		lengths = lengths.astype(int)

		X_Discretization = mydf.copy()
		myModel = hmm.GaussianHMM(n_components=n_components).fit(X, lengths)
		for i in range(0,self.n_gene):
			obs = mydf.iloc[i,:]
			obs.reshape(-1,1)
			Z2 = myModel.predict(obs.reshape(-1,1))
			for j in range(0,self.n_time):
				X_Discretization.iloc[i,j] = Z2[j]

		self.HMM_Z=X_Discretization

	def f_GoI_State(self):
		return self.HMM_Z.loc[self.GoI,:]

	def f_isGoI_same_Z(self):
		if self.HMM_Z.empty:
			self.f_HMM_Discretization()
		isSameState=True
		X_tmp = (self.HMM_Z.loc[self.GoI,:]).T
		X_mark = X_tmp.iloc[:,0]
		for goi in self.GoI:
			if all(X_mark==X_tmp.loc[:,goi]):
				continue
			else:
				isSameState=False	

		return isSameState	


	def f_GoI_same_Z(self):
		if self.f_isGoI_same_Z:
			X_tmp = (self.HMM_Z.loc[self.GoI,:]).T
			X_mark = X_tmp.iloc[:,0]	
			return X_mark.tolist()

	def f_Z_Assoc(self,	Z_cluster = []):
		try:
			X_Prediction = self.HMM_Z.loc[(self.HMM_Z.t0==Z_cluster[0])&
											(self.HMM_Z.t1==Z_cluster[1])&
											(self.HMM_Z.t2==Z_cluster[2])&
											(self.HMM_Z.t3==Z_cluster[3])&
											(self.HMM_Z.t4==Z_cluster[4])&
											(self.HMM_Z.t5==Z_cluster[5])&
											(self.HMM_Z.t6==Z_cluster[6]),:]
		except:
			Z_cluster = self.f_GoI_same_Z()
			X_Prediction = self.HMM_Z.loc[(self.HMM_Z.t0==Z_cluster[0])&
								(self.HMM_Z.t1==Z_cluster[1])&
								(self.HMM_Z.t2==Z_cluster[2])&
								(self.HMM_Z.t3==Z_cluster[3])&
								(self.HMM_Z.t4==Z_cluster[4])&
								(self.HMM_Z.t5==Z_cluster[5])&
								(self.HMM_Z.t6==Z_cluster[6]),:]

		return X_Prediction, X_Prediction.index


	def f_HMM_Prediction(self):
		if self.HMM_Z.empty:
			print 'Call f_HMM_Discretization first to estimate state sequence'

		X_Prediction, glist = self.f_Z_Assoc()

		# Validation 28 genes
		TF_Val = self.TF_Val

		TP_set = list(set(TF_Val)&set(glist))
		Pred_set = list(set(X_Prediction.index) - set(self.GoI))
		# Intersection For Phisical TFs And Prediction
		return TP_set, Pred_set












