from __future__ import print_function

import itertools
import pandas as pd
import numpy as np
import os
import logging
import networkx as nx
import math
from tqdm import tqdm

def buildNetFile(sampdirs, netfile, cutoff, auxDir, writecomponents=False):
    logger = logging.getLogger("grouper")

    sep = os.path.sep
    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampdirs]
    eqfiles = [sep.join([sd, auxDir, '/eq_classes.txt']) for sd in sampdirs]

    tnames = []
    weightDict = {}
    diagCounts = None
    sumCounts = None
    ambigCounts = None
    firstSamp = True
    numSamp = 0
    tot = 0
    eqClasses = {}
    for sffile, eqfile in itertools.izip(sffiles, eqfiles):
        quant = pd.read_table(sffile)
        quant.set_index('Name', inplace=True)

        with open(eqfile) as ifile:
            numSamp += 1
            numTran = int(ifile.readline().rstrip())
            numEq = int(ifile.readline().rstrip())
            logging.info("quant file: {}; eq file: {}; # tran = {}; # eq = {}".format(sffile, eqfile, numTran, numEq))
            if firstSamp:
                for i in xrange(numTran):
                    tnames.append(ifile.readline().rstrip())
                diagCounts = np.zeros(len(tnames))
                sumCounts = np.zeros(len(tnames))
                ambigCounts = np.zeros(len(tnames))
            else:
                for i in xrange(numTran):
                    ifile.readline()

            # for easy access to quantities of interest
            tpm = quant.loc[tnames, 'TPM'].values
            estCount = quant.loc[tnames, 'NumReads'].values
            efflens = quant.loc[tnames, 'EffectiveLength'].values
            epsilon =  np.finfo(float).eps
            sumCounts = np.maximum(sumCounts, estCount)

            for i in xrange(numEq):
                toks = map(int, ifile.readline().rstrip().split('\t'))
                nt = toks[0]
                tids = tuple(toks[1:-1])
                count = toks[-1]
                if tids in eqClasses:
                    eqClasses[tids] += count
                else:
                    eqClasses[tids] = count

                # Add the contribution to the graph
                denom = sum([tpm[t] for t in tids])
                for t1, t2 in itertools.combinations(tids,2):
                    w = count
                    key = (t1, t2)
                    if key in weightDict:
                        weightDict[key] += w
                    else:
                        weightDict[key] = w
                for t in tids:
                    diagCounts[t] += count * (tpm[t] / denom)
                    ambigCounts[t] += count
            firstSamp = False

    lens = quant.loc[tnames, 'Length'].values

    maxWeight = 0.0
    minWeight = 0.0
    prior = 0.1
    edgesToRemove = []

    ##
    #  Go through the weightMap and remove any edges that
    #  have endpoints with too few mapping reads
    ##
    for k,v in weightDict.iteritems():
        c0, c1 = diagCounts[k[0]], diagCounts[k[1]]
        a0, a1 = ambigCounts[k[0]], ambigCounts[k[1]]
        if a0 + a1 > epsilon and a0 > cutoff and a1 > cutoff:
            w = (v+prior) / min((a0+prior), (a1+prior))
            if w > minWeight:
                weightDict[k] = w
                if w > maxWeight:
                    maxWeight = w
            else:
                edgesToRemove.append(k)
        else:
            edgesToRemove.append(k)

    # Actually delete those edges
    for e in edgesToRemove:
        del weightDict[e]

    def nearEnd(tup):
        txp = tup[0]
        pos = tup[1]
        moverhang = 10
        ml = 100
        if pos < -moverhang or pos > lens[txp] + moverhang:
            return False
        elif pos <= ml or pos >= lens[txp] - ml:
            return True
        else:
            return False

    tnamesFilt = []
    relabel = {}
    for i in xrange(len(estCount)):
        if (ambigCounts[i] > cutoff):
            relabel[i] = len(tnamesFilt)
            tnamesFilt.append(tnames[i])
            weightDict[(i, i)] = 1.1

    G = nx.Graph() if writecomponents else None
    with open(netfile, 'w') as ofile:
        writeEdgeList(weightDict, tnames, ofile, G)

    if G is not None:
        clustFile = netfile.split('.net')[0] + '.clust'
        print("Writing connected components as clusters to {}".format(clustFile))
        with open(clustFile, 'w') as ofile:
            cc = nx.connected_component_subgraphs(G)
            for c in cc:
                ofile.write('{}\n'.format('\t'.join(c.nodes())))

def writeEdgeList(weightDict, tnames, ofile, G):
    useGraph = G is not None
    for k,v in weightDict.items():
        ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))
        if useGraph:
            G.add_edge(tnames[k[0]], tnames[k[1]])


def writePajek(weightDict, tnames, relabel, ofile):
    with open(netfile, 'w') as ofile:
        ofile.write("*Vertices\t{}\n".format(len(tnamesFilt)))
        for i, n in enumerate(tnamesFilt):
            ofile.write("{}\t\"{}\"\n".format(i, n))
        ofile.write("*Edges\n")
        print("There are {} edges\n".format(len(weightDict)))
        for k,v in weightDict.items():
            ofile.write("{}\t{}\t{}\n".format(relabel[k[0]], relabel[k[1]], v))

class EquivCollection(object):
    def __init__(self):
        self.tnames = []
        self.eqClasses = {}
        self.hasNames = False

    def setNames(self, names):
        self.tnames = names
        self.hasNames = True

    def add(self, tids, count):
        if tids in self.eqClasses:
            self.eqClasses[tids] += count
        else:
            self.eqClasses[tids] = count

def readEqClass(eqfile, eqCollection):
    with open(eqfile) as ifile:
        numTran = int(ifile.readline().rstrip())
        numEq = int(ifile.readline().rstrip())
        print("file: {}; # tran = {}; # eq = {}".format(eqfile, numTran, numEq))
        if not eqCollection.hasNames:
            tnames = []
            for i in range(numTran):
                tnames.append(ifile.readline().rstrip())
            eqCollection.setNames(tnames)
        else:
            for i in range(numTran):
                ifile.readline()

        for i in range(numEq):
            toks = list(map(int, ifile.readline().rstrip().split('\t')))
            nt = toks[0]
            tids = tuple(toks[1:-1])
            count = toks[-1]
            eqCollection.add(tids, count)

def getCountsFromEquiv(eqCollection):
    countDict = {}
    tn = eqCollection.tnames
    for tids, count in eqCollection.eqClasses.items():
        for t in tids:
            if tn[t] in countDict:
                countDict[tn[t]] += count
            else:
                countDict[tn[t]] = count
    # ensure no division by 0
    for t in eqCollection.tnames:
        if t in countDict:
            countDict[t] += 1.0
        else:
            countDict[t] = 1.0
    return countDict

def flattenClusters(infile, outfile):
    with open(outfile, 'w') as ofile:
        with open(infile) as ifile:
            for i,l in enumerate(ifile):
                toks = l.rstrip().split()
                cname = "cluster{}".format(i)
                for t in toks:
                    ofile.write("{}\t{}\n".format(cname, t))

def filterGraph(expDict, netfile, outfile, auxDir, mincut):
    logger = logging.getLogger("grouper")
    # Get just the set of condition names
    conditions = expDict.keys()
    logging.info("conditions = {}".format(conditions))

    eqClasses = {}
    for cond in conditions:
        print(expDict[cond])
        for sampNum, sampPath in expDict[cond].items():
            if cond not in eqClasses:
                eqClasses[cond] = EquivCollection()
            eqPath = os.path.sep.join([sampPath, auxDir, "/eq_classes.txt"])
            readEqClass(eqPath, eqClasses[cond])

    ambigCounts = {cond : getCountsFromEquiv(eqClasses[cond]) for cond in conditions}

    sailfish = {}
    for cond in conditions:
        sailfish[cond] = ambigCounts[cond]

    G = nx.Graph()
    with open(netfile) as f:
        data = pd.read_table(f, header=None)
        for i in range(len(data)):
            G.add_edge(data[0][i], data[1][i], capacity = float(data[2][i]))

    logging.info("Done reading files for filtering.")
    count = 0
    numTrimmed = 0
    numCut = 0
    cutset = set()
    with open(netfile) as f, open(outfile, 'w') as ofile:
        data = pd.read_table(f, header=None)
        #for i in tqdm(range(len(data))):
        for i in range(len(data)):
            count += 1
            print("\r{} edges checked".format(count), end="")
            #Alternative hypo
            x = data[0][i]
            y = data[1][i]
            non_null=0
            x_all=0
            y_all=0
            for cond in conditions:
                y_g = sailfish[cond][y]
                x_g = sailfish[cond][x]
                r = y_g / x_g
                non_null += (y_g * math.log(r*x_g)) - (r*x_g)
                non_null += (x_g * math.log(x_g)) - x_g
                x_all += x_g
                y_all += y_g
            #null hypothesis
            null = 0
            r_all = y_all / x_all
            for cond in conditions:
                y_g = sailfish[cond][y]
                x_g = sailfish[cond][x]
                mean_x = (x_g + y_g) / (1+r_all)
                null += (y_g * math.log(r_all * mean_x)) - (r_all * mean_x)
                null += (x_g * math.log(mean_x)) - mean_x
            D = 2*(non_null-null)

            if mincut:
                if D > 20 and x != y:
                    s = G.subgraph(nx.shortest_path(G,x))
                    value, partition =  nx.minimum_cut(s, x, y)
                    if value < 10:
                        reachable, non_reachable = partition
                        for e in s.edges_iter(data='capacity'):
                            if (e[0] in reachable and e[1] in non_reachable) or (e[0] in non_reachable and e[1] in reachable):
                                cutset.add((e[0], e[1]))
                        numCut += 1
                    numTrimmed += 1
            else:
                if D > 20 and x != y:
                    numTrimmed += 1
                    cutset.add((x, y))

        G.remove_edges_from(list(cutset))
        for e in G.edges(data = "capacity"):
            ofile.write(e[0] + "\t" + e[1] + "\t" + str(e[2]) + "\n")

    logging.info("Trimmed {} edges".format(numTrimmed))
    logging.info("Cut performed on {} edges".format(numCut))

def addOrphanLinks(sampdirs, auxDir, orphanFileName, cutoff, netFileIn, netFileOut):
    logger = logging.getLogger("grouper")

    sep = os.path.sep
    sffiles = [sep.join([sd, 'quant.sf']) for sd in sampdirs]
    eqfiles = [sep.join([sd, auxDir, '/eq_classes.txt']) for sd in sampdirs]

    tnames = []
    weightDict = {}
    diagCounts = None
    sumCounts = None
    ambigCounts = None
    firstSamp = True
    numSamp = 0
    eqClasses = {}
    for sffile, eqfile in itertools.izip(sffiles, eqfiles):
        quant = pd.read_table(sffile)
        quant.set_index('Name', inplace=True)

        with open(eqfile) as ifile:
            numSamp += 1
            numTran = int(ifile.readline().rstrip())
            numEq = int(ifile.readline().rstrip())
            if firstSamp:
                for i in range(numTran):
                    tnames.append(ifile.readline().rstrip())
                diagCounts = np.zeros(len(tnames))
                sumCounts = np.zeros(len(tnames))
                ambigCounts = np.zeros(len(tnames))
            else:
                for i in range(numTran):
                    ifile.readline()

            # for easy access to quantities of interest
            tpm = quant.loc[tnames, 'TPM'].values
            estCount = quant.loc[tnames, 'NumReads'].values
            efflens = quant.loc[tnames, 'EffectiveLength'].values
            epsilon =  np.finfo(float).eps
            sumCounts = np.maximum(sumCounts, estCount)

            for i in range(numEq):
                toks = map(int, ifile.readline().rstrip().split('\t'))
                nt = toks[0]
                tids = tuple(toks[1:-1])
                count = toks[-1]
                if tids in eqClasses:
                    eqClasses[tids] += count
                else:
                    eqClasses[tids] = count

                denom = sum([tpm[t] for t in tids])
                for t in tids:
                    diagCounts[t] += count * (tpm[t] / denom)
                    ambigCounts[t] += count

            firstSamp = False

    lens = quant.loc[tnames, 'Length'].values
    logging.info("Done reading files for adding orphans")
    # Considering Orphan reads
    def nearEnd(tup):
        txp = tup[0]
        pos = tup[1]
        moverhang = 10
        ml = 100
        if pos < -moverhang or pos > lens[txp] + moverhang:
            return False
        elif pos <= ml or pos >= lens[txp] - ml:
            return True
        else:
            return False

    vertices = set()
    with open(netFileIn, 'r') as ifile:
        for line in ifile:
            line = line.split()
            vertices.add(line[0])
            vertices.add(line[1])

    count = 0
    seenOrphan = {}
    orphanDict = {}
    orphanLinkFiles = [sep.join([sd, auxDir, orphanFileName]) for sd in sampdirs]
    haveLinkFiles = all(os.path.isfile(f) for f in orphanLinkFiles)
    numOrphanLinks = 0
    if haveLinkFiles:
        for olfile in orphanLinkFiles:
            for l in open(olfile):
                left, right = l.rstrip().split(':')
                lp = [map(int, i.split(',')) for i in left.rstrip('\t').split('\t')]
                rp = [map(int, i.split(',')) for i in right.split('\t')]
                lp = [t for t in filter(nearEnd, lp)]
                rp = [t for t in filter(nearEnd, rp)]
                #if len(lp) == 1 or len(rp) == 1:
                for a, b in itertools.product(lp, rp):
                    ltpm = tpm[a[0]] + 10 ** -11  # Laplacian Smoothing
                    rtpm = tpm[b[0]] + 10 ** -11
                    tpm_ratio = 1 - (abs(ltpm - rtpm) / (ltpm + rtpm))
                    read_dist = lens[a[0]] - a[1] + b[1]
                    if tpm_ratio >= .5: # and read_dist <= 300 and tpm[a[0]] > .5 and tpm[b[0]] > .5:
                        a = a[0]; b = b[0]
                        key = (a, b) if a < b else (b, a)
                        if ambigCounts[a] < cutoff or ambigCounts[b] < cutoff:
                            continue
                        c0, c1 = diagCounts[a], diagCounts[b]
                        a0, a1 = ambigCounts[a], ambigCounts[b]
                        if a not in seenOrphan:
                            seenOrphan[a] = set()
                        if b not in seenOrphan:
                            seenOrphan[b] = set()
                        seenOrphan[a].add(b)
                        seenOrphan[b].add(a)
                        if key not in orphanDict:
                            orphanDict[key] = 1.0 / min(a0, a1)
                        else:
                            orphanDict[key] += 1.0 / min(a0, a1)

    for key, value in orphanDict.iteritems():
        if len(seenOrphan[key[0]]) < 3 and len(seenOrphan[key[1]]) < 3:
            if key not in weightDict:
                numOrphanLinks += 1
                weightDict[key] = 1.0 / min(a0, a1)
            else:
                weightDict[key] += 1.0 / min(a0, a1)
                
    logging.info("Added {} orphan link edges".format(numOrphanLinks))

    with open(netFileIn, 'r') as ifile, open(netFileOut, 'w') as ofile:
        for line in ifile:
            ofile.write(line)
        for k,v in weightDict.items():
            ofile.write("{}\t{}\t{}\n".format(tnames[k[0]], tnames[k[1]], v))

    return weightDict;
