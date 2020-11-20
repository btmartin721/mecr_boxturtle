#!/usr/bin/python

import sys
import os
import getopt
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from kneebow.rotor import Rotor
from matplotlib.backends.backend_pdf import PdfPages
from sortedcontainers import SortedDict
from itertools import permutations

def main():
	params = parseArgs()
	
	pd.options.mode.chained_assignment = None  # default='warn'
	
	#read in mu
	mu = loadPickle(params.mu)
	mu_list=list(mu[0].values())[0]
	print("Read mu from", len(mu_list), "replicates.")
	
	#read in sd
	sd = loadPickle(params.sd)
	sd_list=list(sd[0].values())[0]
	print("Read sd from", len(sd_list), "replicates.")	
	
	#parse mu_sd
	dat=parseVAEpickles(mu_list, sd_list)
	#prettyPrint(dat)
	
	o=params.out + ".pdf"
	with PdfPages(o) as pdf:
		#for each replicate, do DBSCAN clustering
		e=list()
		k=list()
		a=list()
		for rep, group in dat.groupby(dat["rep"]):
			#figure out value for epsilon parameter
			if params.f_eps is not None:
				eps=params.f_eps
			else:
				eps=0
				if params.eps=="auto":
					print("\tOptimizing epsilon...",str(rep))
					eps=optimizeEps(group, rep+1)
				else:
					global_sd=np.mean(group["sd"])
					if params.eps=="1sd":
						eps=global_sd
					elif params.eps=="2sd":
						eps=2*global_sd
					elif params.eps=="3sd":
						eps=3*global_sd
					else:
						print("ERROR: Unknown value for eps <-e>:",params.eps)
						sys.exit(1)
			print("Running DBSCAN on replicate",rep)
			assigns=dbscan_cluster(group, eps, params.min_samples)
			group=addPopLabels(group, assigns)
			e.append(eps)
			k.append(len(set(list(group["pop"]))))
			a.extend(group["pop"])
			if params.plot:
				plotClusteredPoints(group, pdf)
		
		if params.plot:
			plotHistogram(k, "K (Clusters)", pdf)
			plotHistogram(e, "DBSCAN eps (epsilon)", pdf)
			
	
	dat["pop"]=a
	oname=params.out+"_assigns.tsv"
	print("Outputting full results for VAE run to:",oname)
	dat.to_csv(oname, sep="\t", header=True, quoting=None, index=False)
	
	props = getProportions(dat)
	if params.matchK:
		print("Running exhaustive packaging-across-K (this might take a while...)")
		props = getProportionsFromT(matchAcrossK(getT(dat)))
	oname=params.out+"_props.tsv"
	print("Outputting VAE assignment proportions to:",oname)
	props.to_csv(oname, sep="\t", header=True, quoting=None, index=False)
	
	
	if params.output_assigns:
		print("<-s> -- Haven't implemented yet... sorry :(")


def getT(dat):
	d=dict()
	for rep, group in dat.groupby(dat["rep"]):
		d["ind"]=list(group["ind"])
		d[rep]=list(group["pop"])
	df=pd.DataFrame.from_dict(d, orient="columns")
	return(df)	

def getProportionsFromT(t):
	l=t.shape[0]
	d=dict()
	d["ind"]=list()
	d["rep"]=list()
	d["pop"]=list()
	for col in t.columns[1:]:
		d["ind"].extend(t["ind"])
		d["rep"].extend([col]*l)
		d["pop"].extend(t[col])
	df=pd.DataFrame.from_dict(d, orient="columns")
	return(getProportions(df))

def matchAcrossK(t):
	#print(t)
	current_best=getMeanBest(getProportionsFromT(t))
	current_t=t
	for rep in t.columns[1:]:
		ord_t=[int(s) for s in set(t[rep])]
		if len(ord_t)==1:
			continue
		perm_bests=list()
		perm_tables=list()
		perm_ords = permutations(ord_t, len(ord_t))
		for p in perm_ords:
			new_col=list(t[rep])
			for value in (ord_t):
				new_col = [str("r"+x) if x==str(value) else x for x in new_col]
			new_col = [str(v[1:]) for v in new_col]
			temp=current_t
			temp[rep]=new_col
			perm_bests.append(getMeanBest(getProportionsFromT(temp)))
			perm_tables.append(temp)
		best=perm_bests.index(max(perm_bests))
		if perm_bests[best] > current_best:
			current_best = perm_bests[best]
			current_t = perm_tables[best]
			print("found better")
	return(current_t)

def getMeanBest(props):
	max=props[props.columns[0:]].max(axis=1)
	return(np.mean(max))

def getProportions(dat):
	df_plot = dat.groupby(['pop', 'ind']).size().reset_index().pivot(columns='pop', index='ind', values=0)
	df_plot=df_plot.fillna(value=0.0)
	return(df_plot)

def dbscan_cluster(group, eps, min_samples):

	points=group[["ae1", "ae2"]].to_numpy()
	#print(points)
	#run DBSCAN
	db=DBSCAN(eps=eps, min_samples=min_samples, algorithm='auto', metric='euclidean').fit(points)
	
	#save cluster labels
	cluster_labels=db.labels_
	#print(cluster_labels)
	
	#build SortedDict to return
	popmap=SortedDict()
	i=0
	for k in group["ind"]:
		pop=str(cluster_labels[i])
		if pop not in popmap:
			l = [k]
			popmap[pop] = l
		else:
			popmap[pop].append(k)
		i+=1
	return(popmap)

def addPopLabels(group, popmap):
	pmap=flattenPopmap(popmap)
	p=list()
	for ind in group["ind"]:
		p.append(pmap[ind])
	group["pop"]=p
	return(group)
	
#function plots clustered coordinates given a SortedDict of coords and a population map
def plotClusteredPoints(group, pdf):
	#set output file name
	plt.figure()
	sns.set(style="ticks")
	#print(group)
	ax = sns.scatterplot(x="ae1", y="ae2", hue="pop", palette="colorblind",data=group)
	#plt.show()	
	plt.plot()
	pdf.savefig()
	plt.close()

#plots a simple histogram
def plotHistogram(dat, name, pdf):
	#of=str(out) + ".snapDistances.pdf"
	sns.set(style="ticks")
	plt.figure()
	x = pd.Series(dat, name=name)
	dp = sns.displot(x, kde=True, rug=True)
	#plt.savefig(of)
	pdf.savefig()
	plt.clf()

#utility function, converts popmap of form key=pop; value=list(inds) to key=ind; value=pop
def flattenPopmap(popmap):
	new_popmap=dict()
	for k in popmap:
		for i in popmap[k]:
			new_popmap[i]=k
	return(new_popmap)

def optimizeEps(group, rep, fig=None):
	"""
	Special thanks to Cory Maklin
	https://towardsdatascience.com/machine-learning-clustering-dbscan-determine-the-optimal-value-for-epsilon-eps-python-example-3100091cfbc
	
	and also stackexchange user georg-un for the kneebow package:
	https://datascience.stackexchange.com/questions/57122/in-elbow-curve-how-to-find-the-point-from-where-the-curve-starts-to-rise
	"""
	X = group[["ae1", "ae2"]].to_numpy()
	neigh = NearestNeighbors(n_neighbors=2)
	nbrs = neigh.fit(X)
	dist, idx = nbrs.kneighbors(X)
	
	dist = np.sort(dist, axis=0)
	d = dist[:,1]
	dist[:,0] = idx[:,0]
	#print(dist)
	#if fig is not None:
	#ax=fig.add_subplot(10,10,rep)
	#ax.plot(d)
	#plt.show()
	
	rotor = Rotor()
	rotor.fit_rotate(dist)
	elbow_index = rotor.get_elbow_index()
	#ax.axhline(dist[elbow_index][1])
	return(dist[elbow_index][1])
	

def scatterPlotReps(dat, out):
	sns.set_theme(style="whitegrid")

	g = sns.relplot(
		data=dat,
		x="ae1", y="ae2",
		hue="ind", size="sd",
		palette="colorblind", sizes=(10, 200),
	)
	g.set(xscale="log", yscale="log")
	g.ax.xaxis.grid(True, "minor", linewidth=.25)
	g.ax.yaxis.grid(True, "minor", linewidth=.25)
	g.despine(left=True, bottom=True)
	#plt.show()

def prettyPrint(df):
	with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
		print(df.sort_values(df.columns[0], ascending=True))

def parseVAEpickles(mu_list, sd_list):
	dat=dict()
	dat["rep"]=list()
	dat["ind"]=list()
	dat["ae1"]=list()
	dat["ae2"]=list()
	dat["sd"]=list()
	for rep in mu_list.keys():
		for idx, coords in enumerate(mu_list[rep]):
			dat["rep"].append(rep)
			dat["ind"].append(idx)
			dat["ae1"].append(coords[0])
			dat["ae2"].append(coords[1])
			dat["sd"].append(np.mean(sd_list[rep][idx]))
	df=pd.DataFrame.from_dict(dat, orient="columns")
	return(df)
	
def loadPickle(pkl):
	ret=list()
	with open(pkl, "rb") as pfh:
		while True:
			try:
				ret.append(pickle.load(pfh))
			except EOFError:
				break
	return(ret)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hm:s:e:i:o:pM:f:sk', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.mu=None
		self.sd=None
		self.eps="2sd"
		self.f_eps=None
		self.out="out"
		self.in_order=None
		self.plot=False
		self.min_samples=1
		self.output_assigns=False
		self.matchK=False


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "h" or opt == "help":
				continue
			elif opt=="m":
				self.mu=arg
			elif opt=="s":
				self.sd=arg
			elif opt=="p":
				self.plot=True
			elif opt=="e":
				if arg not in ["1sd", "2sd", "3sd", "auto"]:
					self.display_help("Invalid value for -e", arg)
				else:
					self.eps=arg
			elif opt=="f":
				self.f_eps=float(arg)
			elif opt=="i":
				self.ind_order=arg
			elif opt=="M":
				self.min_samples=int(arg)
			elif opt=="o":
				self.out=arg
			elif opt=="s":
				self.output_assigns=True
			elif opt=="k":
				self.matchK=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.mu or not self.sd:
			self.display_help("mu and sd files not provided")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nvae_dbscan.py\n")
		print("Author: Tyler Chafin & Bradley Martin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description:Performs DBSCAN clustering given VAE mu and sd along 2 axes")
		print("""
	Mandatory arguments:
		-m	: mu file (.pkl)
		-s	: sd file (.pkl)
		
	Optional arguments:
		-e	: Epislon for DBSCAN
			  Options:
			  1sd     - One global standard deviation (-e 1sd)
			  2sd     - 2*sd (-e 2sd) [default]
			  3sd     - 3*sd (-e 3sd)
			  auto    - Automatically detect optimum eps value
		-f	: Force manually set epsilon (provide value; e.g. -f 0.7)
		-M	: Minimum samples for a DBSCAN cluster [default=1]
		-p	: Output plots using seaborn
		-k	: Exhaustive match-across-K (can take a long time)
		-s	: Generate individual assignment text files per replicate
		-o	: Output file prefix (default=out)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
