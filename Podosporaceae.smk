# -*- snakemake -*-

# from pytools.persistent_dict import PersistentDict
import glob
from Bio.Nexus import Nexus # To read nexus files
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree

### Podosporaceae: Making a phylogeny of the taxa related to Podospora anserina
#############################################################################

# In this study we compiled sequence data for a number of molecular markers
# (ITS, LSU, rpb2, and beta-tubulin) to resolve the phylogenetic relationships
# of lineages within the Podosporaceae family, sensu Wang et al. (2019)
# Studies in Mycology 93:155-252. 
# The primary objective was to determine the relative position of the model
# species *Podospora anserina* and the type species of the genus, *Podospora
# fimiseda*. Unfortunately, there is a lot of missing data.

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/11/26-2020/06/04
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 2

# -------------------------------------------------
## Input data
masternex = config["masternex"]
allmarkersname = config["allmarkersname"]
orders = config["orders"]

# Filtering
minfrac = float(config["minfrac"])
minlen = int(config["minlen"])

# Scripts
sitewise_analyzer = config["sitewise_analyzer"]
genewise_analyzer = config["genewise_analyzer"]
plotShen = config["plotShen"]

## Master lists
outgroup = config["outgroup"]
testmonophyly = config["testmonophyly"]

# -------------------------------------------------

rule all:
	input:
		# This temporal file I don't really need but it's useful to use the checkpoints to produce individual topologies
		"tmp/singlemarkers.tre", 

		# Shen metrics
		"results/ShenMetrics.pdf",

checkpoint nexus2fastas:
	""" Read Nexus into individual markers """
	input:
		aln = masternex
	output:
		directory("singlemarkers")
	run:
		shell("mkdir -p singlemarkers") # The "output"

		# Read nexus file with Bio.Nexus
		master = Nexus.Nexus() # Create a nexus object
		master.read(input.aln) # update

		markers = list(master.charsets.keys()) # Get list of partitions
		print(markers)
		for marker in markers:
			# Get positions included in this partition
			partitioncols = master.charsets[marker]

			# Start an empty msa object
			newalign = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))

			# For each taxon, get those columns:
			for taxon in master.taxlabels:
				record = master.matrix[taxon] # Get sequence of this taxon

				newseq = "" # Start an empty new sequence
				for column in partitioncols:
					newseq += record[column]
				newalign.add_sequence(taxon, newseq) # Add the sliced sequenced to the new msa

			# Write a fasta file for each marker in the working directory
			output_handle = open("singlemarkers/" + marker + ".fas", "w")
			SeqIO.write(newalign, output_handle, "fasta")

# ----------- Single markers
rule minlen:
	""" Remove sequences with too much missing data """
	input:
		fasta = "singlemarkers/{marker}.fas",
	output:
		"cleanmarkers/{marker}" + f".min{minfrac}-{minlen}.fa"
	params:
		threads = 1,
	run:
		# Make a new file
		output_handle = open(output[0], "w")

		for seq_record in SeqIO.parse(input.fasta, "fasta"):
			total = len(seq_record)
			nogaps = len(seq_record.seq.ungap("-"))
			fraction = nogaps/total
			if (fraction >= minfrac) or (nogaps >= minlen):
				SeqIO.write(seq_record, output_handle, "fasta")

		output_handle.close()


rule IQTree_marker:
	""" Make a tree of the marker """
	input:
		alignment = "cleanmarkers/{marker}" + f".min{minfrac}-{minlen}.fa"
	output:
		"iqtree_markers/{marker}" + f".min{minfrac}-{minlen}.treefile",
		"iqtree_markers/{marker}" + f".min{minfrac}-{minlen}.boottrees",
	params:
		bootstraps = 100,
		threads = 4
	shell:
		"iqtree -s {input.alignment} -m MFP -seed 1234 -b {params.bootstraps} -nt {params.threads} -pre iqtree_markers/{wildcards.marker}.min{minfrac}-{minlen}"
		# -b <#replicates>     Bootstrap + ML tree + consensus tree (>=100)
		# -bb Ultrafast bootstraps (UFBoots)
		# -bnni to reduce the risk of overestimating branch supports with UFBoot due to severe model violations.

# --- Put results in a convenient place ---

rule trees_to_results:
	""" Copy the single marker trees into a results folder for convenience """
	input:
		"iqtree_markers/{marker}" + f".min{minfrac}-{minlen}.treefile",
	output:
		"results/{marker}" + f".min{minfrac}-{minlen}.tre",
	shell:
		"cp {input} {output}"

def aggregate_results(wildcards):
	# From Dima's tutorial
	# https://evodify.com/snakemake-checkpoint-tutorial/
	'''
	Aggregate the file names of the random number of files
	generated at the paralogs step
	'''
	checkpoint_output = checkpoints.nexus2fastas.get(**wildcards).output[0] # The folder's name
	return expand("results/{marker}" + f".min{minfrac}-{minlen}.tre",
		marker=glob_wildcards(os.path.join(checkpoint_output, '{marker}.fas')).marker)

rule aggregate_trees:
	""" Collect all trees in a single file """ 
	# With this rule I trigger the formation of all other checkpoint files
	input:
		aggregate_results,
	output:
		"tmp/singlemarkers.tre" # I actually don't need this, but it allows me to get all the trees
	shell:
		"cat {input} > {output}"

## ------ Analysis of Shen et al (2017) statistics ------

def roottree(tree, outgroup, f = 0):
	t = Tree(tree, format = f)
	if len(outgroup) == 1:
		t.set_outgroup(outgroup[0])
	else:
		ancestor = t.get_common_ancestor(outgroup)
		t.set_outgroup(ancestor)
	return(t)

def list2newick(clade): # I didn't use it in the end
	cladestring = '(' + ', '.join(clade) + ')'
	return cladestring

rule make_constrain:
	""" Produce contrasting topologies to evaluate the T1 = (A,B) vs T2 = (A,C) dispute """
	input:
		tree = "results/{marker}" + f".min{minfrac}-{minlen}.tre"
	output:
		constrain = "ShenMetrics/{marker}" + f".min{minfrac}-{minlen}-constrain.tre"
	run:
		print(input.tree)

		# What outgroup taxa survived?
		t = Tree(input.tree)
		leaves = [node.name for node in t.iter_leaves()]  # same as [leaf for leaf in t.iter_leaf_names()]
		
		# Make a list with only the outgroups that survived the filtering
		survoutgroup = [taxon for taxon in outgroup if taxon in leaves]

		# Find what tree is this, T1 or T2?
		t = roottree(input.tree, survoutgroup) # Root first
		t.prune(testmonophyly) # remove all taxa that are not relevant for the test
		print(t)

		# Write the tree in the topology that already has
		if t.check_monophyly(values=[testmonophyly[0], testmonophyly[1]], target_attr="name")[0]: # AB
			constrain = Tree(f"(({testmonophyly[0]}, {testmonophyly[2]})), ({testmonophyly[1]}), {list2newick(survoutgroup)};", format=1) # constrain tree is AC

		elif t.check_monophyly(values=[testmonophyly[0], testmonophyly[2]], target_attr="name")[0]: # AC
			constrain = Tree(f"(({testmonophyly[0]}, {testmonophyly[1]})), ({testmonophyly[2]}), {list2newick(survoutgroup)};", format=1) # constrain tree is AB

		elif t.check_monophyly(values=[testmonophyly[1], testmonophyly[2]], target_attr="name")[0]: # BC:
			constrain = Tree(f"(({testmonophyly[0]}, {testmonophyly[1]})), ({testmonophyly[2]}), {list2newick(survoutgroup)};", format=1) # constrain tree is AB

		constrain.write(format = 9, outfile = output.constrain)


rule IQTREE_constrained_single:
	""" Make an alternative tree with the topological constrain """
	input:
		alignment = "cleanmarkers/{marker}" + f".min{minfrac}-{minlen}.fa",
		tree = "ShenMetrics/{marker}" + f".min{minfrac}-{minlen}-constrain.tre"
	output:
		"ShenMetrics/{marker}" + f".min{minfrac}-{minlen}-constrained.treefile"
	params:
		bootstraps = 100,
		threads = 3
	shell:
		"iqtree -s {input.alignment} -m MFP -nt {params.threads} -g {input.tree} -pre ShenMetrics/{wildcards.marker}.min{minfrac}-{minlen}-constrained"


rule RAxML_SLS:
	""" Use RAxML to calculate the site-wise log-likelihood score """
	input:
		alignment = "cleanmarkers/{marker}" + f".min{minfrac}-{minlen}.fa",
		T1 = "results/{marker}" + f".min{minfrac}-{minlen}.tre",
		T2 = "ShenMetrics/{marker}" + f".min{minfrac}-{minlen}-constrained.treefile"
	output:
		trees = "ShenMetrics/T1_T2_{marker}.tre",
		sls = "ShenMetrics/RAxML_perSiteLLs.{marker}_podospora_site_lk"
	shell:
		"cat {input.T1} {input.T2} > {output.trees}; "
		"cd ShenMetrics; "
		"raxmlHPC -f G -s ../{input.alignment} -m GTRGAMMA -z T1_T2_{wildcards.marker}.tre -n {wildcards.marker}_podospora_site_lk"


def get_main_alignment(wildcards):
	""" Get name of relevant alignment for the Shen et al analysis """
	checkpoint_output = checkpoints.nexus2fastas.get(**wildcards).output[0] # I don't really need this but it activates the checkpoints
	return expand("ShenMetrics/RAxML_perSiteLLs.{marker}_podospora_site_lk",
		marker=allmarkersname)

rule Shen2017:
	""" Run the scripts from Shen et al 2017 """
	input:
		sls = get_main_alignment,
		orders = orders,
	output:
		"ShenMetrics/podospora_sitwise_statistics.txt",
		"ShenMetrics/podospora_genewise_lnL.txt",
	shell:
		"cp {input.orders} ShenMetrics/; "
		"cd ShenMetrics; "
		"perl ../{sitewise_analyzer} && perl ../{genewise_analyzer}" 

rule plotShen2017:
	""" Plot SLS and dGLS """
	input: 
		sitwise = "ShenMetrics/podospora_sitwise_statistics.txt",
		genewise = "ShenMetrics/podospora_genewise_lnL.txt",
	output:
		"results/ShenMetrics.pdf"
	script:
		plotShen
