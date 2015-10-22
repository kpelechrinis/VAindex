#!/usr/bin/env python

import sys
import numpy as np
import random
import os

mu_l = [-21, 21];
v_l = [0.2, 0.6, 1.5, 5];
cov_l = [0, 0.25, 0.6, 1];
no_types = 16
dim = 5
no_instances = 1500

for delta in [0, 0.2, 0.4, 0.6, 0.8, 1]:
	os.system("mkdir data/delta_"+str(delta))
	for svar in [0, 1, 2]:
		os.system("mkdir data/delta_"+str(delta)+"/svar_"+str(svar))
		for cor in [0, 1, 2]:
			print delta,svar,cor
			os.system("mkdir data/delta_"+str(delta)+"/svar_"+str(svar)+"/cor_"+str(cor))
			f_prototype = "data/delta_"+str(delta)+"/svar_"+str(svar)+"/cor_"+str(cor)+"/prototypes.txt" 
                        f_pr = open(f_prototype,"w") 
			print >> f_pr,"TypeID","\t","Mean Vector","\t","Cov Matrix"
			f_nodes = "data/delta_"+str(delta)+"/svar_"+str(svar)+"/cor_"+str(cor)+"/nodes.txt"
			f_n = open(f_nodes,"w")
			print >> f_n,"nodeID","\t","Attribute Vector"
			for types in range(no_types):
				x_mean = [0, 0, 0, 0, 0];
				for i in range(len(x_mean)):
					x_mean[i] = mu_l[0] + (mu_l[1]-mu_l[0])*random.random()
				cov = [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]
				for i in range(len(x_mean)):
					cov[i][i] = pow(v_l[svar]*abs(x_mean[i])+(v_l[svar+1]-v_l[svar])*abs(x_mean[i])*random.random(),2)
				for i in range(len(x_mean)):
					for j in range(i+1,len(x_mean)):
						if random.random() < delta:
							#c = np.sign(-1+2*random.random())
                                                        c = 1
							cov_f = cov_l[cor]+(cov_l[cor+1]-cov_l[cor])*random.random()
                                                        cov[i][j] = np.sqrt(cov[i][i])*np.sqrt(cov[j][j])*c*cov_f
						else:
							cov[i][j] = 0
						cov[j][i] = cov[i][j]
				u, v = np.linalg.eig(cov)
				while (min(u) < 0):
					print "Not a positive semi-definite covariance matrix!"
					for i in range(len(x_mean)):
                                        	cov[i][i] = pow(v_l[svar]*abs(x_mean[i])+(v_l[svar+1]-v_l[svar])*abs(x_mean[i])*random.random(),2)
                                	for i in range(len(x_mean)):
                                        	for j in range(i+1,len(x_mean)):
                                                	if random.random() < delta:
                                                        	#c = np.sign(-1+2*random.random())
                                                        	c = 1
                                                        	cov_f = cov_l[cor]+(cov_l[cor+1]-cov_l[cor])*random.random()
                                                        	cov[i][j] = np.sqrt(cov[i][i])*np.sqrt(cov[j][j])*c*cov_f
                                                	else:
                                                        	cov[i][j] = 0
                                                	cov[j][i] = cov[i][j]
					u, v = np.linalg.eig(cov)
				print >> f_pr,types,"\t",x_mean,"\t",cov
				for instance in range(no_instances):
					x_tmp = np.random.multivariate_normal(x_mean,cov)
					print >> f_n,str(instance)+"-"+str(types),"\t",str(x_tmp[0])+","+str(x_tmp[1])+","+str(x_tmp[2])+","+str(x_tmp[3])+","+str(x_tmp[4])	
			f_pr.close()
			f_n.close()
