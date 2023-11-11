import numpy as np
import pandas as pd
from functools import partial
from collections import defaultdict
from sklearn import metrics
import src.load_data as ld
import src.set_analysis_func as func


if __name__=='__main__':
	parser = argparse.ArgumentParser(
		descripition='calculate gene sets similarity in a embedding space')
	parser.add_argument('--emb', dest='emb_f', type=str,
						help='input file path and file name for embedding')
	parser.add_argument('--genelist', dest='genelist_f', type=str,
						help='input file path and file name for embedding genes')
	parser.add_argument('--geneset', dest='geneset_f', type=str,
						help='input file path and file name for the gene set database')
	parser.add_argument('--rankedlist', dest='rankedlist_f', type=str,
						help='input file path and file name for the ranked gene list')
	parser.add_argument('--out', dest='out_f', type=str,
						help='output file path and name')
	parser.add_argument('-n', '--n_processor', dest='n_process', type=int, default=10
						help='number of processors')
	parser.add_argument('--min', dest='min_size', type=int, default=10
						help='the minimum number of genes in a set')
	parser.add_argument('--max', dest='max_size', type=int, default=300
						help='the maximum number of genes in a set')

	args = parser.parse_args()

	# load embedding
	node_vectors = np.loadtxt(emb_f, delimiter=',')
	node_list = []
	with open(genelist_f, 'r') as f:
	    for line in f:
	        node_list.append(line.strip())
	        
	S = metrics.pairwise.cosine_similarity(node_vectors, node_vectors)

	# create gene to embedding id mapping
	g_node2index = {j:i for i,j in enumerate(node_list)}
	g_index2node = {i:j for i,j in enumerate(node_list)}
	g_node2index = defaultdict(lambda:-1, g_node2index)


	# load gene set data
	geneset = ld.load_gmt(geneset_f)

	geneset_indices = ld.term2indexes(
		geneset, g_node2index, upper=max_size, lower=min_size)


	# get background genes
	geneset_all_genes = set()
	for x in geneset:
	    set_all_genes = set_all_genes.union(geneset[x])
	    
	geneset_all_genes = geneset_all_genes.intersection(node_list)
	geneset_all_indices = [g_node2index[x] for x in geneset_all_genes]
	geneset_terms = list(geneset_all_indices.keys())




	#load ranked list
	ranked_list = pd.read_csv(rankedlist_f, sep='\t', 
		index_col=0, header=None)
	ranked_list = [str(y) for y in ranked_list.index]
	ranked_list = [g_node2index[y] for y in ranked_list if y in node_list]


	f = partial(func.gsea_andes, ranked_list=ranked_list, matrix=S, 
	            term2indices=geneset_indices, annotated_indices=geneset_all_indices)

	with Pool(n_process) as p:
		rets = p.map(f, geneset_terms)

	zscores = pd.DataFrame(rets, index=geneset_terms)
