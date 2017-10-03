import sys
import os
import argparse
import itertools
import statistics
import re
import csv
import fileinput
import numpy as np
from collections import defaultdict
from subprocess import call
import pandas as pd
from collections import Counter

alpha = 0.8 #ratio alloted to the probability value for labels
newEdgeProb = 0.9 #add a new edge if the shared label prob is greater than this
                    
def getMedWeight(graph, node1, node2):
	weights = []
	for (x, weight) in graph[node1]:
		if weight != 1.1:
			weights.append(weight)
	for (x, weight) in graph[node2]:
		if weight != 1.1:
			weights.append(weight)
	
	if not weights:
		return(0)
	else:
		return(statistics.median(weights))

def readLabels(keys, iterNum):
	#below we get labelToContigs = labels -> contigs that map to this label
	#and contigToLabels = contigs -> all labels
	#contigLabelsToProb = contig, labels -> prob
	labelToContigs = defaultdict(set)
	contigToLabels = defaultdict(set)
	contigLabelsToProb = {}
	with open(keys["seed_file"], 'r') as ifile:
		for line in ifile:
			data = (line.strip('\n')).split('\t')
			contigName = data[0]
			if (iterNum >= 1):
				data = data[3].split()
				curLabel = ""
				for i in range(1, len(data), 2):
					curLabel = data[i-1]
					if (curLabel != "" and curLabel != "__DUMMY__" and data[i].lower() != "nan" and curLabel != "IGNORE"):
						contigToLabels[contigName].add(curLabel)
						labelToContigs[curLabel].add(contigName)
						contigLabelsToProb[(contigName, curLabel)] = float(data[i])
			else:
				curLabel = data[1]
				contigToLabels[contigName].add(curLabel)
				labelToContigs[curLabel].add(contigName)
				contigLabelsToProb[(contigName, curLabel)] = float(data[i])
	return(labelToContigs, contigToLabels, contigLabelsToProb)

#get the sum of the probabilities that 2 nodes have the same label
#@profile
def getProb(contigToLabels, contigLabelsToProb, node1, node2):
	totalProbability = 0
	sharedLabels = contigToLabels[node1].intersection(contigToLabels[node2])
	for label in sharedLabels:
		totalProbability += (contigLabelsToProb[(node1, label)] * contigLabelsToProb[(node2, label)])
	return totalProbability

#@profile
def changeEdgeWeights(orgGraph, graph, contigToLabels, contigLabelsToProb, ofile):
	changesMade = 0
	noChange = 0
	weightCalc = 0
	selfLoops = 0

	for (node1, node2), weight in graph.iteritems():
		if node1 == node2:
			ofile.write(node1 + "\t" + node2 + "\t1.1\n")
			selfLoops += 1
		else:
			prob = getProb(contigToLabels, contigLabelsToProb, node1, node2)
			if prob > 0:
				for (x, w) in orgGraph[node1]:
					if x == node2:
						orgWeight = w
				for (x, w) in orgGraph[node2]:
					if x == node1:
						orgWeight = w
				newWeight = ((1 - alpha) * orgWeight) + (prob * alpha)
				changesMade += 1
			else:
				newWeight = weight
				noChange += 1
			weightCalc += newWeight
			ofile.write(node1 + "\t" + node2 + "\t" + str(newWeight) + "\n")
	weightCalc /= (changesMade+noChange)

	return (weightCalc, (changesMade + noChange + selfLoops))

#@profile
def addNewEdges(orgGraph, graph, contigToLabels, labelToContigs, contigLabelsToProb, ofile):
	changesMade = 0
	weightCalc = 0
	for label in labelToContigs:
		for node1, node2 in itertools.combinations(labelToContigs[label], 2):
			if node1 != node2:
				if (node1, node2) not in graph and (node2, node1) not in graph:
					totalProb = getProb(contigToLabels, contigLabelsToProb, node1, node2)
					if (totalProb > newEdgeProb):
						orgWeight = getMedWeight(orgGraph, node1, node2)
						newWeight = ((1 - alpha) * orgWeight) + (totalProb * alpha)
						ofile.write(node1 + "\t" + node2 + "\t" + str(newWeight) + "\n")
						orgGraph[node1].add((node2, newWeight))
						changesMade += 1
						weightCalc += newWeight
	if (changesMade > 0):
		weightCalc /= changesMade
	return (weightCalc, changesMade)

def run(keys, labelFile, juntoConfigFile, netFile, outdir):
	
	if (netFile == keys["graph_file"]):
		print ("ERROR: junto graph file should be different from the net file otherwise it will be overwritten")
		return 0

	orgGraph = defaultdict(set)
	orgGraphSize = 0
	avgGraphWeight = 0
	with open(netFile, 'r') as ifile:
		for line in ifile:
			edge = (line.strip('\n')).split('\t')
			orgGraph[edge[0]].add((edge[1], float(edge[2])))
			orgGraphSize += 1
			if (float(edge[2]) != 1.1):
				avgGraphWeight += float(edge[2])
	avgGraphWeight /= orgGraphSize

	diff = 10000
	prevDiff = 0
	iters = 1
	i = 0

	with open(keys["graph_file"], 'w') as f:
		for node1 in orgGraph:
			for node2, weight in orgGraph[node1]:
				if node1 == node2:
					f.write(node1 + "\t" + node2 + "\t" + str(0.00001) + "\n")
				else:
					f.write(node1 + "\t" + node2 + "\t" + str(weight) + "\n")

	#while (abs(prevDiff - diff) > (0.1*prevDiff)):
	while (diff > (0.05*prevDiff)):
		i += 1
		print ("Started iteration number: " + str(i) + "\n")

		if (i == 2):
			prevDiff = diff;

		tempfile = open("temp.txt", 'w')
		with open(keys["graph_file"], 'r') as ifile:
			for line in ifile:
				edge = (line.strip('\n')).split('\t')
				node1 = edge[0]
				node2 = edge[1]
				weight = float(edge[2])
				if node1 == node2:
					tempfile.write(node1 + "\t" + node2 + "\t" + str(0.00001) + "\n")
				else:
					tempfile.write(node1 + "\t" + node2 + "\t" + str(weight) + "\n")
			
		tempfile.close()
		call(["mv", "temp.txt", keys["graph_file"]])

		call(["junto", "config", juntoConfigFile])
		call(["cp", keys["output_file"], keys["seed_file"]])

		print("Done running label propogation: processing results.")
		graph = {}
		graphSize = 0
		with open(keys["graph_file"], 'r') as ifile:
			for line in ifile:
				edge = (line.strip('\n')).split('\t')
				if edge[0] != edge[1]:
					graph[(edge[0], edge[1])] = float(edge[2])
				else:
					graph[(edge[0], edge[1])] = 1.1
				graphSize += 1

		labelToContigs, contigToLabels, contigLabelsToProb = readLabels(keys, i)
		print("Done reading labels.")
		sizeNewGraph = 0

		with open(keys["graph_file"], 'w') as ofile:
			#write previous edges with new weights
			(avgOldWeight, temp) = changeEdgeWeights(orgGraph, graph, contigToLabels, contigLabelsToProb, ofile)
			sizeNewGraph += temp
			#add new edges between contigs with same labels
			(avgNewWeight, temp) = addNewEdges(orgGraph, graph, contigToLabels, labelToContigs, contigLabelsToProb, ofile)
			sizeNewGraph += temp

		oldGraphWeight = avgGraphWeight
		avgGraphWeight = (avgOldWeight + avgNewWeight)/2
		
		with open(keys["seed_file"], 'w') as ofile:
			for contig in contigToLabels:
				curProb = 0
				for label in contigToLabels[contig]:
					temp = contigLabelsToProb[(contig, label)]
					ofile.write(contig + "\t" + label + "\t" + str(temp) + "\n")
					curProb += temp
				if (curProb < 1):
					ofile.write(contig + "\t" + "IGNORE" + "\t" + str(1-curProb) + "\n")

		diff = abs(sizeNewGraph - graphSize)
		print("Size diff: {}".format(diff))

		with open(juntoConfigFile, "w") as configFile:
                	for x in keys:
                        	configFile.write(x + " = " + keys[x] + "\n")
                	configFile.write("data_format = edge_factored\n")
                	configFile.write("iters = " + str(iters) + "\n")
               		configFile.write("prune_threshold = 0\n")
                	configFile.write("algo = adsorption\n")
