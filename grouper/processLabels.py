import os
from subprocess import call
import pandas as pd
from collections import defaultdict
import logging

#run BLAST and return list of label files in order [TS->ASdb, AS->TSdb] 
def runBLAST(faFile1, faFile2, outdir, threads):
	logger = logging.getLogger("grouper")
	TSdb = os.path.sep.join([outdir, "TS", "db"])
	ASdb = os.path.sep.join([outdir, "AS", "db"])
	call(["makeblastdb", "-in", faFile1, "-dbtype", "nucl", "-out", TSdb])
	call(["makeblastdb", "-in", faFile2, "-dbtype", "nucl", "-out", ASdb])

	labelFile1 = os.path.sep.join([outdir, "TS.ASdb" ]) 
	labelFile2 = os.path.sep.join([outdir, "AS.TSdb" ]) 
	#run blastn for TS against AS
	call(["blastn", "-db", ASdb, "-query", faFile1, "-outfmt", "6", "-out", labelFile1, "-num_threads", str(threads)])
	call(["blastn", "-db", TSdb, "-query", faFile2, "-outfmt", "6", "-out", labelFile2, "-num_threads", str(threads)])
	
	logging.info("Done running BLAST using {} threads".format(threads))
	return [labelFile1, labelFile2]

#get the lengths of transcripts from one of the quant files
def getLengths(expDict):
    logger = logging.getLogger("grouper")
    lengths = {}
    condition = list(expDict)[0]
    quantFile = os.path.sep.join([expDict[condition][0], "quant.sf"])
    logging.info("Getting txp lengths from: {}".format(quantFile))

    lengths = pd.read_table(quantFile, names=["Name", "Length"], usecols=["Name", "Length"])
    lengths.set_index("Name", inplace=True)
    lengths = lengths.to_dict()['Length']

    return lengths

#read in the BLAST results and process them to get two way best match. This creates the seed file for Junto.
#Label files should be in order [TS->ASdb, AS->TSdb] 
def genFinalLabels(expDict, keys, labelFiles, finalLabelFile, outdir):
	lengths = getLengths(expDict)

	blastOut = open(keys["seed_file"], 'w')
	if (len(labelFiles)==1):
		with open(labelFiles[0].strip(), 'r') as f:
			data = pd.read_table(f, header=None, usecols=[0,1], names=['query', 'subject'])
			table1 = data.set_index("query").to_dict()['subject']
		try:
			for key,value in table1.items():
				probability=1.0
				blastOut.write(key + "\t" + value + "\t" +str(probability)+ "\n")
		except KeyError:
			pass
	else:
		labelFiles[0] = labelFiles[0].strip()
		labelFiles[1] = labelFiles[1].strip()
		outfile1 = outdir + "best" + os.path.basename(labelFiles[0])
		outfile2 = outdir + "best" + os.path.basename(labelFiles[1])
		command = "cat " + labelFiles[0]  + " |  sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > " + outfile1
		call(command, shell=True)
		command = "cat " + labelFiles[1] + " |  sort -k2,2 -k12,12nr -k11,11n | sort -u -k2,2 --merge > " + outfile2
		call(command, shell=True)
		with open(outfile1, 'r') as f:
			table1 = {}
			lentable1 = {}
			for line in f:
				line = line.split()
				table1[line[0]] = line[1]
				lentable1[line[0]] = abs(float(line[6]) - float(line[7]))

		with open(outfile2, 'r') as f:
			table2 = {}
			lentable2 = {}
			for line in f:
				line = line.split()
				table2[line[1]] = line[0]
				lentable2[line[1]] = abs(float(line[8]) - float(line[9]))	

		call(["rm", outfile1])
		call(["rm", outfile2])

		try:
			for key,value in table1.items():
				if key in table2:
					alignLen = max(lentable1[key], lentable2[key])
					probability=min(1.0, (alignLen / float(lengths[key])))
					blastOut.write(key + "\t" + value + "\t" +str(probability)+ "\n")
		except KeyError:
				pass
	blastOut.close()
	call(["cp", keys["seed_file"], finalLabelFile])
