#### Boas Pucker ###
### Luna Theresa Lingemann ###
### pucker@uni-bonn.de ###
### v0.4###

__usage__ = """
					python heatmap_plotter.py
					--genes <GENES_FILE>
					--exp <EXPRESSION_FILE>
					--out <FIGURE_OUTPUT_FILE>
					
					optional:
					--norm (activates normalization per gene)
					"""

import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
# --- end of imports --- #


def load_genes( gene_file ):
	"""! @brief load genes """
	
	genes = {}
	gene_order = []
	with open( gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip()
			if "\t" in parts:
				parts = parts.split('\t')
				genes.update( { parts[0]: parts[1] } )
				gene_order.append( parts[0] )
			else:
				genes.update( { parts: parts } )
				gene_order.append( parts )
			line = f.readline()
	return genes, gene_order


def load_exp( gene_file ):
	"""! @brief load expression values """
	
	exp = {}
	with open( gene_file, "r" ) as f:
		sample_order = f.readline().strip().split('\t')[1:]
		for header in sample_order:
			exp.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			values = map( float, parts[1:] )
			for idx, val in enumerate( values ):
				exp[ sample_order[ idx ] ].update( { parts[0]: val } )
			line = f.readline()
	return exp, sample_order


def generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_order, norm ):
	"""! @brief generate figure """
	
	exp_dataset = []	#[ [ column1 ], [column2], [column3] ]
	gene_lables = []
	for gene in gene_order:
		values = []
		for sample in sample_order:
			values.append( exp[ sample ][ gene ] )
		exp_dataset.append( values )
		gene_lables.append( genes[ gene ] )
	data = pd.DataFrame( exp_dataset, columns=sample_order )
	data.index = gene_lables + []
	
	if norm:
		data_z = data.apply( lambda x: (x-x.mean())/x.std(), axis = 1)
	else:
		data_z = data

	fig, ax = plt.subplots(figsize=(20, 16))
	
	h = sns.heatmap( data_z, annot=data, fmt=".1f", cmap=LinearSegmentedColormap.from_list("custom_heatmap", ["#efefef", "#d5a6bd", "#c27ba0"]), annot_kws={"size": 10, "weight": "bold"}, xticklabels=1, yticklabels=1, cbar_kws={"label": "Z-Score"})

	new_names = ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"]
	h.set_xticklabels( new_names, fontsize=12, fontweight='bold', fontstyle='italic')
	ax.set(xlabel="", ylabel="")
	ax.xaxis.tick_top()
	h.set_yticklabels( h.get_yticklabels(), rotation=0, fontweight='bold', fontstyle='italic')	#turn y-axis labels to horizontal
	cbar = h.collections[0].colorbar
	plt.setp(cbar.ax.get_yticklabels(), fontsize=14, fontweight="bold")
	cbar.ax.set_ylabel("Z-Score", fontsize=20, fontweight="bold")
	
	fig.tight_layout()
	
	#plt.subplots_adjust( left=0.2, right=0.999, top=0.9, bottom=0.01, hspace=0.03  )
	
	fig.savefig( figfile, dpi=300 )
	plt.close( "all" )


def main( arguments ):
	"""! @brief run generation of plots """
	
	gene_file = arguments[ arguments.index('--genes')+1 ]
	exp_file = arguments[ arguments.index('--exp')+1 ]
	figfile = arguments[ arguments.index('--out')+1 ]
	
	if '--norm' in arguments:
		norm = True
	else:
		norm = False
	
	genes, gene_order = load_genes( gene_file )

	exp, sample_order = load_exp( exp_file )
	if '--samples' in arguments:
		sample_order = arguments[ arguments.index('--samples')+1 ].split(',')

	generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_order, norm )


if '--genes' in sys.argv and '--exp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
