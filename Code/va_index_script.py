#!/usr/bin/env python

import sys
import numpy as np
import random
import math
import os
import igraph
from sklearn.metrics.pairwise import cosine_similarity
from scipy import spatial
from scipy.spatial.distance import pdist

def my_euclidean_similarity(x,y,gamma):
	temp = np.array(x) - np.array(y)
	return math.exp(-np.linalg.norm(temp)/gamma)

def my_cosine_similarity(x,y):
	return (sum(np.array(x)*np.array(y)))/(np.linalg.norm(x)*np.linalg.norm(y))

def my_similarity(x,y,metric):
	return pdist(np.array([x,y]),metric)[0]

p_b = 0.3

f_results = open("assortativity_results_euclidean_pow_tmp.txt","w")

metric = 'euclidean'
power = 1

gamma = ["gamma_0","gamma_0.2","gamma_0.8","gamma_0.6","gamma_0.4","gamma_1"]
svar = ["svar_0","svar_1","svar_2"]
cor = ["cor_0","cor_1","cor_2"]
no_types = [2, 4, 8, 16]

for mcsim in range(5000):
 print "Iteration: ",mcsim
 for g in [random.choice(gamma)]:
    for s in [random.choice(svar)]:
        for c in [random.choice(cor)]:
            # we need to create a dictionary with the vector attributes
            f_vec = open(g+"/"+s+"/"+c+"/"+"nodes.txt","r")
            vecs = dict()
            for line in f_vec:
                line.rstrip()
                fields = line.rsplit("\t")
                vecs[fields[0].rstrip()]= []
                if "-" in fields[0]:
                    elements = fields[1].rsplit(',')
                    for e in range(len(elements)):
                        vecs[fields[0].rstrip()].append(float(elements[e].rstrip()))
	    f_vec.close()
	    for t in [random.choice(no_types)]:
				for p_t in [random.random()]:
                    			for i in range(1):
						#we need to identify how many nodes of each type we will have
						#total number of nodes in the network is 1000
						break_points = sorted(random.sample(range(1000),t-1),key=int)
						node_IDs = list()
						#we need to pick which t of the 16 types we will use
						types = sorted(random.sample(range(16),t),key=int)
						#type_perm will keep permutaitons of the IDs of different types 
						type_perm = list()
						for j in range(16):
							type_perm.append(np.random.permutation(1500))
						#create the node_IDs that will participate in the network
						count = 0
						for j in range(len(break_points)):
							count_type = 0;
							while (count < break_points[j]):
								node_IDs.append(str(type_perm[types[j]][count_type])+"-"+str(types[j]))
								count = count + 1
								count_type = count_type + 1
                        #we need to add the last type
                        			count_type = 0;
                        			while (count < 1000):
                            				node_IDs.append(str(type_perm[types[-1]][count_type])+"-"+str(types[-1]))
                            				count = count + 1
							count_type = count_type + 1
						f_out = open("tmp.txt","w")
						for node1 in range(1000):
							for node2 in range(node1+1,1000):
								#check if they are of the same type
								if (node_IDs[node1].rsplit('-')[1] == node_IDs[node2].rsplit('-')[1]):
									if (random.random() <= p_t):
										print >> f_out,node_IDs[node1]+" "+node_IDs[node2]
								else:
									if (random.random() <= p_b):
										print >> f_out,node_IDs[node1]+" "+node_IDs[node2]
                        			f_out.close()
                        			# Calculate ground truth assortativity & baseline
                        			graph = igraph.Graph.Read_Ncol("tmp.txt")
                        			types = []
						types_num = dict()
						for dim in ['0','1','2','3','4']:
							types_num[dim] = []
						vector_assort = []
                        			for n in range(len(graph.vs['name'])):
                            				types.append(int(graph.vs[n]['name'].rsplit('-')[1]))
							for dim in ['0','1','2','3','4']:
								types_num[dim].append(vecs[graph.vs[n]['name']][int(dim)])
                        			true_assort = graph.assortativity_nominal(types,directed=False)
						for dim in ['0','1','2','3','4']:
							vector_assort.append(graph.assortativity(types_num[dim],directed=False))
						mean_assort = np.mean(vector_assort)
                        			# Calculate VA-index
                        			# The network structure we created is now in tmp.txt
                        			f_net = open("tmp.txt","r")
                        			edge_count = 0
						similarity_real = []
						distance_real = []
                        			nodes_tmp = dict()
                        			for edge in f_net:
                            				edge_count = edge_count + 1
                            				edge.rstrip()
                            				nodes = edge.rsplit(" ")
							distance_real.append(my_similarity(vecs[nodes[0].rstrip()],vecs[nodes[1].rstrip()],metric))
                            				if nodes[0].rstrip() not in nodes_tmp:
                                    				nodes_tmp[nodes[0].rstrip()] = 1
                            				if nodes[1].rstrip() not in nodes_tmp:
                                    				nodes_tmp[nodes[1].rstrip()] = 1
                        			avg_similarity_rand = []
						sim_list = []
						dis_list = []
						max_dist = 0
						for n1 in nodes_tmp.keys():
							for n2 in nodes_tmp.keys():
								if n1 != n2:
									tmp = my_similarity(vecs[n1],vecs[n2],metric)
									dis_list.append(tmp)
									if (tmp > max_dist) and (tmp>0):
										max_dist = tmp
						if metric == 'braycurtis':
							max_dist = 1
						similarity_real=list(1-np.power(np.array(distance_real)/max_dist,power))
						sim_list = list(1-np.power(np.array(dis_list)/max_dist,power))
						avg_similarity_real = np.mean(similarity_real)
						std_rand = []
						max_d = []
                        			for rand_net in range(50):
                            				tmp_sim = 0
							indices = random.sample(range(len(sim_list)),edge_count)
							tmp_sim = sum(list(np.array(sim_list)[indices]))
							std_rand.append(np.std(list(np.array(sim_list)[indices])))
							max_d.append(np.max(list(np.array(dis_list)[indices])))
                            				avg_similarity_rand.append(tmp_sim/edge_count)
						pooled_var = np.mean(std_rand)
                        			if (avg_similarity_real <= np.percentile(avg_similarity_rand,2.75)) or (avg_similarity_real >= np.percentile(avg_similarity_rand,97.5)):
                            				empirical_p = "(p-value < 0.05)"
                        			else:
                            				empirical_p = "(p-value > 0.05)"
                        			f_net.close()
						print >> f_results, g+" "+s+" "+c+" "+str(t)+" "+str(p_t)+"\t",true_assort,"\t",avg_similarity_real - np.mean(avg_similarity_rand),"\t",pooled_var,"\t",mean_assort
						os.system("rm tmp.txt")

f_results.close() 
