# Grouper 

Grouper is a tool for clustering and annotating contigs from *de novo* transcriptome assemblies.  There are two main modules in Grouper: the clustering module and the labeling module.  The former is based on the tool, [RapClust](https://github.com/COMBINE-lab/RapClust), and is designed to be run downstream of the [Sailfish](https://github.com/kingsfordgroup/sailfish) or [Salmon](https://github.com/COMBINE-lab/salmon) tools for rapid transcript-level quantification.  It relies on the *fragment equivalence classes*, orphaned read mappings and quantification information computed by these tools in order to determine how contigs in the assembly are potentially related and cluster them accordingly.  The labeling module in Grouper is able to incorporate information from annotated genomes of closely related species to annotate contigs in the *de novo* assembly.  Hence, the different modules of Grouper are able to efficiently utilize information from multiple sources to accurately cluster and annotate contigs from transcriptome assemblies. 

## Dependencies


The clustering module of Grouper depends on the [MCL](http://micans.org/mcl/) clustering tool (to be available in the environment where it runs).

Similarly, the labeling module depends on the [Junto](https://github.com/parthatalukdar/junto) library for label propagation (to be available in the environment where it runs). This will require the relevant Java version. You can add this by cloning the repository and running the following commands:

```
export JUNTO_DIR=<path to junto folder>
export PATH="$PATH:$JUNTO_DIR/bin"
```

Further, Grouper depends on the following Python packages:
  
  1. [Click](http://click.pocoo.org/5/)
  2. [PyYAML](https://pypi.python.org/pypi/PyYAML)
  3. [Pandas](http://pandas.pydata.org/)
  4. [NumPy](http://www.numpy.org/)
  5. [Networkx v1.11](https://networkx.github.io/)

However, you should be able to install Grouper via `pip` and have these python dependencies installed automatically.  To install Grouper via pip, you can use:

```
> pip install biogrouper
```

You should now have a `Grouper` executable in your path.  You can test this with the following command:

```
> Grouper --help
```

You should see the following output:

```
Usage: Grouper [OPTIONS]

Options:
  --config TEXT  Config file describing the experimental setup
  --help         Show this message and exit.
```

### Notes on Grouper options

Grouper provides flags to enable certain features that may improve the accuracy of clustering.  Specifically, one can use information from orphaned reads to link contigs prior to clustering (setting `orphan` to true in the YAML file), and can apply a min cut filter to force the clustering module to separate certain contigs based on their expression statistics across conditions.  It is not the case, however, that applying these features always improves the accuracy of Grouper.  While there is no single simple rule that one can follow to decide if one or both of these options should be applied, we observed the following when testing Grouper on a wide range of data.  The more complete and accurate the assembly, the less useful these features tend to be.  Hence, if your assembly is of high quality (as evalutated by e.g. [Transrate](http://hibberdlab.com/transrate/) or [Detonate](http://deweylab.biostat.wisc.edu/detonate/)), then you may wish to run Grouper without these flags enabled, as it may yield superior performance.  However, if your assembly is of lower quality, and contains many orphaned reads, then you may wish to enable these flags.  We are exploring a more thorough methodology for determining which set of options is optimal for a given data set, and will update this documentation as we determine more concrete recommendations. 

## Using Grouper

Grouper is written in Python and is easy to use.  Below, we explain how to use it with Salmon.  There are two main steps involved in running Grouper:

  1. Run Salmon on each sample in your experiment, passing it the `--dumpEq` option.  This will tell Salmon to dump a representation of the fragment equivalence classes that it computed during quasi-mapping of each sample.  If you wish to use orphan read information for joining contigs in Grouper, use the `--writeOrphanLinks` option as well, which will dump orphan read pair information to a file.  Apart from these additional option, Salmon should be run normally (i.e. passing in whatever other options are appropriate for your samples).  
  2. Run Grouper, providing it with a configuration file that describes the experimental setup of your samples, and where the Salmon quantification results have been written. You can also choose whether or not to use the additional filters and the labeling module in Grouper.
    
Let's illustrate this pipeline with a particular example, the following experimental data from the [Trapnell et al. paper](http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html):

Accession | Condition | Replicate
----------|-----------|----------
SRR493366 | scramble  | 1
SRR493367	| scramble  | 2
SRR493368	| scramble  | 3
SRR493369	| HOXA1KD	  | 1
SRR493370	| HOXA1KD	  | 2
SRR493371 | HOXA1KD   | 3

We'll assume that the raw read files reside in the directory `reads`.  Assuming that you've already built the index on the transcriptome you wish to quantify, a typical run of Salmon on this data would look something like.

```
> parallel -j 6 "samp={}; salmon quant -i index -l a -1 <(gunzip -c reads/{$samp}_1.fq.gz) -2 <(gunzip -c reads/{$samp}_2.fq.gz) -o {$samp}_quant --dumpEq --writeOrphanLinks -p 4" ::: SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
```

This will quantify each sample, and write the result to the directory `samplename_quant`.  Given this setup, we're now ready to run the clustering module in Grouper.  First, we have to make an appropriate config file.  We demonstrate one using both the optinal filters in Grouper:

```
conditions:
    - Control
    - HOXA1 Knockdown
samples:
    Control:
        - SRR493366_quant
        - SRR493367_quant
        - SRR493368_quant
    HOXA1 Knockdown:
        - SRR493369_quant
        - SRR493370_quant
        - SRR493371_quant
outdir: human_grouper
orphan: True
mincut: True
```

you can place this in a file called `config.yaml`.  Grouper uses [YAML](http://yaml.org/) to specify its configuration files.  The configuration file must contain the following three entries; `conditions`, `samples`, and `outdir`.  The `conditions` entry lists the conditions present in the sample. The `samples` entry is a nested dictionary of lists; there is a key corrseponding to each condition listed in the `conditions` entry, and the value associated with this key is a list of quantification directories of the samples for this condition.  Finally, the `outdir` entry specifies where the Grouper output and intermediate files should be stored.  Optionally, the `orphan` and `mincut` entries tell Grouper which extra filters to use. If these lines are not added to the config file, by default, the filters are not applied. Given the above, we can run Grouper as:

```
> Grouper --config config.yaml
```

This will process the samples, generate the mapping ambiguity graph, filter it according to the conditions and the optional filter, and cluster the resuling graph (Grouper uses [MCL](http://micans.org/mcl/) internally for clustering).  Once Grouper is finished, the `human_grouper` directory should exist.  It will contain the following files:

`mag.clust, mag.filt.net,  mag.flat.clust,  mag.net,  stats.json, log.txt, mag.orphan.net`

The most important file for downstream processing is `mag.flat.clust`.  It contains the computed cluster information in a "transcript-to-gene" mapping formation that is compatible with downstream tools like [tximport](https://github.com/mikelove/tximport).  The other files may be useful for exploration, but they are more intended for Grouper's internal use (e.g. `mag.filt.net` contains the filtered mapping ambiguity graph that is used for clustering).

## Labeling Module
-------------------

In order to annotate contigs in the assembly, the labeling module of Grouper requires information from a closely related species.  Our test species in this example is human and the closely related annotated species is chimp.  This information can be added to the config file in in *any* **one** of the following formats:

  1. You can pass the FASTA files to Grouper in the following way and it will run a two-way BLAST assigning seed labels to contigs. Ensure that the FASTA files are passed in the following order (the first is from the test species, second from the annotated species)

  ```
  fasta:
      - human.transcripts.fa
      - chimp.transcripts.fa
  ```

  2. If you have already run BLAST, you can pass the output files (in BLAST outfmt 6). Again, ensure that the first one is BLAST of contigs from test species against the annotated species and the second is BLAST of contigs from annotated species against the test species. 

  ```
  labels:
      - human.chimpdb.txt
      - chimp.humandb.txt
  ```

  3. If you wish to use a pre-processed label file, you can pass a two-column file where the first is the set of contigs from the test species and second the label. If a contig has multiple labels in the input file, one will be chosen arbitrarily as seed.

  ```
  labels:
      - human.labels.txt
  ```

So a sample config file provided with the FASTA files (example 1) would look something like this:
```
conditions:
    - Control
    - HOXA1 Knockdown
samples:
    Control:
        - SRR493366_quant
        - SRR493367_quant
        - SRR493368_quant
    HOXA1 Knockdown:
        - SRR493369_quant
        - SRR493370_quant
        - SRR493371_quant
fasta:
    - human.transcripts.fa
    - chimp.transcripts.fa
outdir: human_grouper
orphan: True
mincut: True
threads: 12
```

This also uses the optional filters in Grouper to generate the mapping ambiguity graph and runs BLAST using 12 threads (if this is not specified, it is run using 8 threads by default). The ouput directory in this case will contain a sub-folder `Annotated` with the following files:

`final.labels.txt, label.graph.txt, label.mag.clust, label.mag.flat.clust, label.stats.json, raw.label.graph.txt, seed.labels.txt`

Once again, the `mag.flat.clust` contains the computed cluster information in a "transcript-to-gene" mapping formation that can be used for downstream analyses. The file `seed.labels.txt` contains the initial contig to gene labeling. More importantly, the file `final.labels.txt` contains the labels after running Grouper and a contig may have multiple labels in this file, each with an associated score. The rest of the files are for internal use in the algorithm.

## Citations:
-------------

Experiments in Graph-based Semi-Supervised Learning Methods for Class-Instance Acquisition. Partha Pratim Talukdar, Fernando Pereira, ACL 2010

Differential analysis of gene regulation at transcript resolution with RNA-seq by Cole Trapnell, David G Henderickson, Martin Savageau, Loyal Goff, John L Rinn and Lior Pachter, Nature Biotechnology 31, 46–53 (2013).

Stijn van Dongen. Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht, 2000

Charlotte Soneson, Michael I Love, and Mark D Robinson. Differential analyses for rna-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4, 2015.
