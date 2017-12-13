# t-SNE [1] is a tool to visualize high-dimensional data. 
# It converts similarities between data points to joint probabilities 
# and tries to minimize the Kullback-Leibler divergence between the joint probabilities
 # of the low-dimensional embedding and the high-dimensional data. 
 # t-SNE has a cost function that is not convex, 
 # i.e. with different initializations we can get different results.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from ggplot import *
from sklearn.manifold import TSNE
import itertools
import pandas as pd
from Operations import Operations


class Visualization(object):
	def __init__(self, df=pd.DataFrame(), GoI=[], Dir = 'Img'):
		if not df.empty:
			self.df = df
		else:
			db = pd.read_csv('myDB_Validation.csv',index_col = 'id')
			X1 = db.loc[:,'t0':'t6']
			X2 = db.loc[:,'TF':'Chen_isTar']
			GoI = ['AT2G37090','AT4G18780','AT5G15630','AT5G17420','AT5G44030'] 
			X_GoI = X1.loc[GoI,'t0':'t6']
			X_pTFs = X1.loc[X2.TF==1,:]
			X = pd.concat([X_GoI, X_pTFs])
			self.df = X

		self.GoI = GoI
		self.df_tsne=[]
		self.ngene = len(df)
		self.label = []
		self.df_label = [];
		self.Dir = Dir
		
		if not (os.path.isdir(Dir)):
			os.mkdir(Dir)

	def f_GoI(self,List_GoI):
		self.GoI = List_GoI

	def f_setlabel(self, label):
		self.label = label

	def f_set_df_tsne(self, df_tsne):
		self.df_tsne = df_tsne

	def f_set_df_label(self, df_label):
		if len(df_label) == self.ngene:
			self.df_label = df_label
			# self.label = df_label[:,-1]
			self.label = df_label.label
		else:
			print 'Caution: Dimension Of Label Incompatible v.s. Input Dataframe!'
			raw_input('Press To Continue...')

	def f_normalization(self, df,isTime=True,isObj=False):
		if isTime:
			dfs = df.T
			dfs = ( dfs - dfs.mean() ) / dfs.std()
			return dfs

	def f_plot_trajectory(self, isLegend=True, isGoI=False,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Trajectories.pdf'):
		dfs = self.f_normalization(self.df, isTime = True)
		# print dfs
		gh = dfs.plot(
			title='Trajectories'
			, lw=3 , legend=isLegend, style=style )
		if isLabel:
			# plot background
			dfs = self.f_normalization(self.df, isTime = True)
			gh = dfs.plot(
				title='Trajectories'
				, lw=3 , legend=False, style = 'lightgrey')

			gh = dfs.plot(
				title='Trajectories'
				, lw=3 , legend=False, style = 'lightgrey')
		if isGoI:
			# plot background
			dfs = self.f_normalization(self.df, isTime = True)
			gh = dfs.plot(
				title='Trajectories'
				, lw=3 , legend=False, style = 'lightgrey')
			
			# plot frontground
			try:
				df_goi = self.df.loc[self.GoI,:]		
				df_goi = self.f_normalization(df_goi, isTime = True)
				df_goi.plot(lw=4,ax=gh, legend=isLegend)
			except:
				print 'No GoI temperal profiles are available.'

		gh.set_xlim(xmin=0)
		gh.set_xlim(xmax=6)
		gh.set_xlabel("Time (hour)")
		gh.set_ylabel("Normalized Expression Level")
		gh.set_xticklabels(['0','0.25','0.5','1','4','12','24'])
		# gh.legend(loc = 'best')

		# Output Figures
		fig = plt.gcf()
		plt.grid(isGrid)
		fig.set_size_inches(18.5,10.5)
		if isOutput:
			fig.savefig(self.Dir+'/'+FileName,bbox_inches='tight')

	def f_get_TSNE(self):
		df = self.f_normalization(self.df, isTime = True)
		df = df.T
		tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=5000,random_state=0)
		tsne_results = tsne.fit_transform(df)
		self.df_tsne = df.copy()
		self.df_tsne['x-tsne'] = tsne_results[:,0]
		self.df_tsne['y-tsne'] = tsne_results[:,1]

	def f_plot_TSNE(self, FileName='TSEN.pdf', isTrans=True, isLabel=False, isGoI=False, isOutput=True):
		df_color = self.df_tsne.copy()
		df_color['label'] = 'DEGs'

		if isTrans:
			chart = ggplot( df_color, aes(x='x-tsne', y='y-tsne'))\
					 + geom_point(size=70,alpha=0.8) \
					 + ggtitle("tSNE dimensions")
			chart



		if isLabel:
			df_color['label'] = self.label
			df_color['label'] = df_color['label'].apply(lambda i: 'Class'+str(i))

			chart = ggplot( df_color, aes(x='x-tsne', y='y-tsne', color='label'))\
				 + geom_point(size=70,alpha=0.8) \
				 + ggtitle("tSNE dimensions") #\
				 # + scale_color_brewer(type = 'qual')
			chart


		if isGoI:
			df_color = self.df_tsne.copy()
			df_color['label'] = 'Genes'
			for goi in self.GoI:
				df_color.loc[goi,'label'] = goi
			df_color = df_color.sort_values(by='label',ascending=False)

			chart = ggplot( df_color, aes(x='x-tsne', y='y-tsne', color='label'))\
					 + geom_point(size=70,alpha=0.8) \
					 + ggtitle("tSNE dimensions") \
					 # + scale_color_brewer(type = 'qual')
			chart	
		if isOutput:
			chart.save(self.Dir+'/'+FileName)


if __name__ == '__main__':
	print('----------Plot Signatures----------')
	ops = Operations()
	ops.f_getSignature()

	vis = Visualization(df=ops.sig_timecourse)
	vis.f_plot_trajectory(isLegend=False, isGoI=False,
		isOutput=True,isGrid=False, isLabel=False, style=None, FileName='Sig_timecourse.pdf')

	vis.f_set_df_tsne(ops.sig_tSNE)
	vis.f_plot_TSNE(FileName='Sig_tSNE.pdf', isTrans=True, isLabel=False, isGoI=False, isOutput=True)














