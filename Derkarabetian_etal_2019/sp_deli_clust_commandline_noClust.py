#!/usr/bin/env python
# coding: utf-8

#########################################################################################
#########################################################################################
## Script from:
## Derkarabetian S., Castillo S., Peter K.K., Ovchinnikov S., Hedin M.
## 2019.
## "An Empirical Demonstration of Unsupervised Machine Learning in Species Delimitation".
## Molecular Phylogenetics and Evolution, 139: 106562.
## https://doi.org/10.1016/j.ympev.2019.106562.
#########################################################################################
#########################################################################################

#########################################################################################
## Script has been modified from the original.
## We implemented an early stopping callback to mitigate VAE overfitting,
## and we also modified it to perform a user-specified number of replicate
## runs for each datasets with various missing data and minor allele frequency filters.
## The output from the replicate runs is written to pickle objects and PDFs.
#########################################################################################

# Load needed libraries

import argparse
import pickle
import sys

import keras
from keras.layers import *
from keras.models import Sequential,Model
from keras import backend as K
from keras.callbacks import EarlyStopping

import tensorflow as tf

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

def Get_Arguments():
	"""
	Parse command-line arguments. Imported with argparse.
	Returns: object of command-line arguments.
	"""
	parser = argparse.ArgumentParser(description="Run VAE species delimitation for varied datasets", add_help=False)

	required_args = parser.add_argument_group("Required Arguments")
	optional_args = parser.add_argument_group("Optional Arguments")
	
	## Required Arguments
	required_args.add_argument("-r", "--runs",
								type=int,
								required=True,
								help="Number of runs per dataset.")

	required_args.add_argument("-i", "--ind",
								type=str,
								required=True,
								help="Max missing data percent for individual filtering datasets")

	required_args.add_argument("-p", "--pop",
								type=str,
								required=True,
								help="Max missing data percent for population)filtering datasets")

	required_args.add_argument("-m", "--maf",
								type=str,
								required=True,
								help="Minor allele frequency setting to run.")
	
	optional_args.add_argument("-e", "--epochs",
								type=int,
								required=False,
								default=1000,
								help="Max number of epochs to run. Default: 1000")

	optional_args.add_argument("-s", "--samples",
								type=int,
								required=False,
								default=100,
								help="Number of random samples used to make std. deviation blobs.")
	
	optional_args.add_argument("-h", "--help",
								action="help",
								help="Displays this help menu")

	if len(sys.argv)==1:
		print("\nExiting because no command-line options were called.\n")
		parser.print_help(sys.stderr)
		sys.exit(1)

	args = parser.parse_args()
	return args


# Functions to save and load datasets

def save_data(filename,tmp):
	fo = open(filename, "w")
	for n in range(tmp["shape"][0]):
		print(tmp["name"][n],tmp["group"][n],end="",file=fo)
		for i in range(tmp["shape"][1]):
			print("",end=" ",file=fo)
			print(",".join([str(i) for i in list(tmp["one_hot"][n][i])]),end="",file=fo)
		print("",file=fo)
	fo.close()
	
def load_data(filename, test_size=0.3):
	tmp = {
		"data_set":filename.replace(".txt",""),
		"name":[],
		"group":[],
		"one_hot":[]
  	}
	file = open(filename,"r")
	for line in file:
		val = line.rstrip("\n").split(" ")
		tmp["name"].append(val[0])
		tmp["group"].append(val[1])
		tmp["one_hot"].append([])
	
		for n in range(2,len(val)):
			tmp["one_hot"][-1].append(val[n].split(","))
	  
	tmp["name"] = np.array(tmp["name"],np.str)
	tmp["group"] = np.array(tmp["group"],np.str)
	tmp["one_hot"] = np.array(tmp["one_hot"],np.float)
	tmp["train"], tmp["val"] = train_test_split(tmp["one_hot"],test_size=test_size) # split into training, test datasets
	tmp["shape"] = tmp["one_hot"].shape
	tmp["train_shape"] = tmp["train"].shape
	tmp["val_shape"] = tmp["val"].shape

	return tmp


# Functions for plot generation

def plot_z(z,loc,colors,title="none",legend=True):
	x = z[:,0]
	y = z[:,1]
	colorcount = 0  

	for lo in np.unique(loc):
		xx = []
		yy = []
		for i in range(z.shape[0]):
			if lo == loc[i]:
				xx.append(x[i])
				yy.append(y[i])
		plt.scatter(xx,yy,label=lo,c=colors[colorcount])
		colorcount += 1

	if title != "none": plt.title(title)
	if legend == True:
		plt.legend(bbox_to_anchor=(1, 0, 0.5, 1),loc="upper left",)
  
def plot_mu_sg(mu,sg,loc,colors,sample=100,alpha=0.5,title="none",legend=True):
	colorcount = 0
	for lo in np.unique(loc):
		xx = []
		yy = []
		for i in range(len(loc)):
			if lo == loc[i]:
				x,y = (mu[i] + np.random.normal(0,1,size=(sample,2)) * sg[i]).T
				xx += list(x)
				yy += list(y)
		plt.scatter(xx,yy,label=lo,alpha=alpha,c=colors[colorcount])
		colorcount += 1

	plt.scatter(mu[:,0],mu[:,1],s=30,facecolors='none',edgecolors='black',linewidths=1)
	if title != "none": plt.title(title)
	if legend == True:
		plt.legend(bbox_to_anchor=(1, 0, 0.5, 1),loc="upper left",)


# Function to make VAE model

def mk_model(
	original_dim,    # number of snps
	cat,             # number of categories of one-hot-encoding
	latent_dim=2,
	
	# encoder
	en_dim=[100,100,100],
	en_drop=[0.5,0.5,0.5],
	
	# decoder
	de_dim=[100,100,100],  # number of neurons for each layer
	de_drop=[0.5,0.5,0.5], # rate of dropout for each layer
	
	act = "elu" # activation function for each layer
):
  
	def act_fn(fn,tensor):
		if fn == "leakyrelu": return LeakyReLU()(tensor)
		else: return Activation(fn)(tensor)
  
	half_cat = int(cat/2)
	########################################################################
	# INPUT
	########################################################################  
	x_in = Input(shape=(original_dim,cat),name="x_in")
	x_in_em = Dense(half_cat,use_bias=False,name="x_in_em")(x_in)
	
	########################################################################
	# ENCODER :: Q(z|X)
	########################################################################
	en = Flatten()(x_in_em)
	en = BatchNormalization(scale=False,center=False)(en)
  
	for i in range(len(en_dim)):
		en = Dense(en_dim[i])(en)
		en = Dropout(en_drop[i])(en)
		en = act_fn(act,en)
		en = BatchNormalization(scale=False,center=False)(en)
	
	########################################################################
	# Z (Latent space)
	########################################################################  
	Z_mu = Dense(latent_dim)(en)
	Z_log_sigma_sq = Dense(latent_dim)(en)
  
	Z_sigma = Lambda(lambda x: K.exp(0.5*x))(Z_log_sigma_sq)
  
	Z = Lambda(lambda x: x[0]+x[1]*K.random_normal(K.shape(x[0])))([Z_mu,Z_sigma])
  
	########################################################################
	# DECODER :: P(X|z)
	########################################################################
	de = Z
	for i in range(len(de_dim)):
		de = Dense(de_dim[i])(de)
		de = Dropout(de_drop[i])(de)
		de = act_fn(act,de)
		de = BatchNormalization(scale=False,center=False)(de)
  
	de = Dense(original_dim*half_cat)(de)
	########################################################################
	# OUTPUT
	########################################################################  
	x_out_em = Reshape((-1,half_cat))(de)
	x_out = Dense(cat,activation="softmax")(x_out_em)
	########################################################################
  
	def vae_loss(kl_weight=0.5):
		def loss(x_true, x_pred):
			# mask out missing data!
			mask = K.sum(x_in,axis=-1)

			# sigma (or standard deviation), keeping close to 1
			# mu (or mean), keeping close to 0
			kl_loss = K.sum(K.square(Z_mu) + K.square(Z_sigma) - Z_log_sigma_sq - 1.0, axis=-1)

			# reconstruction (categorical crossentropy)
			recon = K.sum(keras.metrics.categorical_crossentropy(x_in,x_out) * mask, axis=-1)
	  
			return K.mean(recon + kl_loss * kl_weight)
		return loss
	
	def acc(x_true,x_pred):
		mask = K.sum(x_in,axis=-1,keepdims=True)
		acc = K.sum(K.square(x_in-x_out),axis=-1,keepdims=True)
		return K.mean(1.0 - K.sqrt(K.sum(acc*mask,axis=1)/K.sum(mask,axis=1)))
	
	vae0 = Model([x_in],[x_out],name="vae0")
	vae0.compile(optimizer='adam', loss=vae_loss(0.1), metrics=[acc])
  
	vae1 = Model([x_in],[x_out],name="vae1")
	vae1.compile(optimizer='adam', loss=vae_loss(0.5), metrics=[acc])
  
	enc = Model([x_in],[Z_mu,Z_sigma],name="enc")
	return vae0,vae1,enc


# function to generate models, train, and plot the latent space (Z)

def do_it(data, epoch_setting, run_number, pdf):
	plt.rcParams['figure.figsize'] = [10, 10]
	plt.style.use('seaborn-colorblind')

	def gen(dat, batch_size, key):
		while True:
			idx = np.random.randint(0,len(dat[key])-1,size=batch_size)
			tmp = dat[key][idx]
			yield tmp,tmp

	K.clear_session()
	vae0,vae1,enc = mk_model(data["shape"][1],data["shape"][2])
	loss_history = []
	acc_history = []
	val_loss_history = []
	val_acc_history = []
  
	# Stop when loss starts to climb
	earlystopping = EarlyStopping(monitor='val_loss', patience=50, restore_best_weights=True)
  
	# we do different batch_sizes: 1/4, 1/3, 1/2 and 1/1 of the data
	# (similar to changing the temperature)
	r = 4
	for i in range(r):
		f = 1/(r-i)
		batch_size = int(data["shape"][0] * f + 0.5)
		steps = int(data["shape"][0]/batch_size + 0.5)
		epochs = int(epoch_setting * f + 0.5)

		for vae in (vae0,vae1):
			print("-")
			his = vae.fit_generator(
				gen(data, batch_size, "train"),
				steps_per_epoch=steps,
				epochs=epochs,
				verbose=False, 
				validation_data=gen(data, batch_size, "val"), 
				validation_steps=steps,
				callbacks=[earlystopping] # stop when loss starts to climb
			)
			loss_history += list(his.history['loss'])
			acc_history += list(his.history['acc'])
			val_loss_history += list(his.history['val_loss'])
			val_acc_history += list(his.history['val_acc'])
		
	# Set colors for my specific dataset. Requires modification if you have different number of groups.
	colors = ["gray", "cyan", "darkorange", "yellow", "magenta", "blue", "red", "rebeccapurple", "forestgreen"]

	if i == r-1:
		
		fig = plt.figure(run_number)
		plt.subplot(2, 2, 1)
		plt.plot(np.arange(len(loss_history)),loss_history)
		plt.plot(np.arange(len(val_loss_history)), val_loss_history)
		plt.xlabel("Epoch", fontsize=14)
		plt.ylabel("Loss", fontsize=14)
		plt.title(data["data_set"])
		plt.legend(["Training", "Test"], loc="best", prop={'size': 14})
		plt.subplot(2, 2, 2)
		plt.plot(np.arange(len(acc_history)),acc_history)
		plt.plot(np.arange(len(val_acc_history)), val_acc_history)
		plt.legend(["Training", "Test"], loc="best", prop={'size': 14})
		plt.ylabel("Accuracy", fontsize=14)
		plt.xlabel("Epoch", fontsize=14)

		vae_mu,vae_sg = enc.predict(data["one_hot"])
		plt.subplot(2, 2, 3)
		plot_z(vae_mu,data["group"],colors,legend=False)
		plt.subplot(2, 2, 4)
		plot_mu_sg(vae_mu,vae_sg,data["group"],colors,sample=100)
		pdf.savefig(fig)
	  
	return(vae_mu,vae_sg,colors)


def calc_mu_sd(mu,sg,loc,colors,sample=100,alpha=0.5,title="none",legend=True):
	colorcount = 0
	grp_x = list()
	grp_y = list()
	for lo in np.unique(loc):
		xx = []
		yy = []
		for i in range(len(loc)):
			if lo == loc[i]:
				x,y = (mu[i] + np.random.normal(0,1,size=(sample,2)) * sg[i]).T
				xx += list(x)
				yy += list(y)
				grp_x += xx
				grp_y += yy
  
	return grp_x, grp_y
	

def do_sil(X, max_clust=9):
	
	sil = dict()
	for i in range(1, max_clust):
		
		n_clusters = i+1
		
		# Create a subplot with 1 row and 2 columns
		fig, (ax1, ax2) = plt.subplots(1, 2)
		fig.set_size_inches(7, 7)

		# The 1st subplot is the silhouette plot
		# The silhouette coefficient can range from -1, 1 but in this example all
		# lie within [-0.1, 1]
		ax1.set_xlim([-0.1, 1])
		# The (n_clusters+1)*10 is for inserting blank space between silhouette
		# plots of individual clusters, to demarcate them clearly.
		ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

		# Initialize the clusterer with n_clusters value and a random generator
		# seed of 10 for reproducibility.
		clusterer = KMeans(n_clusters=n_clusters)
		cluster_labels = clusterer.fit_predict(X)

		# The silhouette_score gives the average value for all the samples.
		# This gives a perspective into the density and separation of the formed
		# clusters
		silhouette_avg = silhouette_score(X, cluster_labels)
		print("For n_clusters =", n_clusters,
			  "The average silhouette_score is :", silhouette_avg)
		sil[n_clusters] = silhouette_avg
		
		# Compute the silhouette scores for each sample
		sample_silhouette_values = silhouette_samples(X, cluster_labels)

		y_lower = 10
		for i in range(n_clusters):
			# Aggregate the silhouette scores for samples belonging to
			# cluster i, and sort them
			ith_cluster_silhouette_values =                 sample_silhouette_values[cluster_labels == i]

			ith_cluster_silhouette_values.sort()

			size_cluster_i = ith_cluster_silhouette_values.shape[0]
			y_upper = y_lower + size_cluster_i

			color = cm.nipy_spectral(float(i) / n_clusters)
			ax1.fill_betweenx(np.arange(y_lower, y_upper),
							  0, ith_cluster_silhouette_values,
							  facecolor=color, edgecolor=color, alpha=0.7)

			# Label the silhouette plots with their cluster numbers at the middle
			ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

			# Compute the new y_lower for next plot
			y_lower = y_upper + 10  # 10 for the 0 samples

		ax1.set_title("The silhouette plot for the various clusters.")
		ax1.set_xlabel("The silhouette coefficient values")
		ax1.set_ylabel("Cluster label")

		# The vertical line for average silhouette score of all the values
		ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

		ax1.set_yticks([])  # Clear the yaxis labels / ticks
		ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

		# 2nd Plot showing the actual clusters formed
		colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
		ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
					c=colors, edgecolor='k')

		# Labeling the clusters
		centers = clusterer.cluster_centers_
		# Draw white circles at cluster centers
		ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
					c="white", alpha=1, s=200, edgecolor='k')

		for i, c in enumerate(centers):
			ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
						s=50, edgecolor='k')

		ax2.set_title("The visualization of the clustered data.")
		ax2.set_xlabel("Feature space for the 1st feature")
		ax2.set_ylabel("Feature space for the 2nd feature")

		plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
					  "with n_clusters = %d" % n_clusters),
					 fontsize=14, fontweight='bold')
 
	return sil

def do_clust(data, run_number, current_mu, current_sg, samples):
	colors = ["gray", "cyan", "darkorange", "yellow", "magenta", "blue", "red", "rebeccapurple", "forestgreen"]
	grp_x, grp_y = calc_mu_sd(current_mu, current_sg, data["group"], colors, sample=samples)

	# Convert values to numpy arrays and do column_stack.
	x = np.array(grp_x)
	y = np.array(grp_y)
	tmp = np.column_stack((x, y))
	a = np.array(list(tmp))

	sil_avg = do_sil(a)
	
	return sil_avg

#########################################################################
## MAIN ##
#########################################################################

print("tensorflow version: {}".format(tf.__version__))
print("keras version: {}".format(keras.__version__))

args = Get_Arguments()

maf = args.maf

ind = args.ind
pop = args.pop

sd_samples = int(args.samples)

mu_pkl = "mu_missInd{0}_Pop{1}_maf{2}.pkl".format(ind, pop, maf)
sg_pkl = "sg_missInd{0}_Pop{1}_maf{2}.pkl".format(ind, pop, maf)

epoch_setting = int(args.epochs)
runs = int(args.runs)

data_set = "missInd{0}_Pop{1}_maf{2}.onehot.txt".format(ind, pop, maf)

mu_dict = dict()
sg_dict = dict()
mu_tmp = dict()
sg_tmp = dict()
print("\n\n\nWorking on dataset {}".format(data_set))
with PdfPages("{}.pdf".format(data_set)) as pdf:
	for i in range(runs):
		print("\tRun {}".format(i+1))
		mu, sg, colors = do_it(load_data("data/"+data_set, test_size=0.2), epoch_setting, i+1, pdf)
		mu_tmp[i] = mu
		sg_tmp[i] = sg

mu_dict[data_set] = mu_tmp
sg_dict[data_set] = sg_tmp

with open(mu_pkl, "wb") as mufile:
	pickle.dump(mu_dict, mufile)
	
with open(sg_pkl, "wb") as sgfile:
	pickle.dump(sg_dict, sgfile)
	
print("\nDone!\n")
sys.exit(0)

