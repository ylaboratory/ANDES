# ANDES

## Algorithm for Network Data Embedding and Similarity analysis

ANDES is a suite of standalone scripts for comparing similarity between gene sets using precomputed gene embeddings.
It includes a consensus protein–protein interaction network embedding generated with node2vec and a sample
geneset database (Gene Ontology Biological Process genesets for _Homo Sapiens_).

ANDES has two main functions:

  1. calculating gene sets similarity in a embedding space
  2. ranked-based GSEA using gene embedding information

These functions are implemented in `src/set_analysis_fun.py` and the demo
jupyter notebook (`demo.ipynb`) shows sample usage.

### Features

- Gene-set similarity: Compute pairwise similarity scores between two gene sets in embedding space.
- Embedding-based GSEA: Perform a ranked Gene Set Enrichment Analysis (GSEA) using embedding-derived gene rankings.

## Citation

If you use ANDES in your work, please cite:
> [A best-match approach for gene set analyses in embedding spaces.](https://pubmed.ncbi.nlm.nih.gov/39231608/)
Li L, Dannenfelser R, Cruz C, Yao V. Genome Research. 2024.

## Installation

1. Install conda if you haven't already
2. Create and activate the ANDES environment:
   
```sh
conda env create -f env.yml
conda activate ANDES
```

## Usage

To quickly get started we recommend looking at our `demo.ipynb`. Alternatively, ANDES can be run
from the command line in both modes with the following commands.

Compute similarity between all pairs of genesets in two databases / gmt files:

```sh
python src/andes.py --emb embedding_file.csv --genelist embedding_gene_ids.txt --geneset1 first_gene_set_database.gmt --geneset2 second_gene_set_database.gmt --out output_file.csv -n num_processor
```

Compute a ranked-based comparison for a geneset database (such as Gene Ontology) given a ranked list of genes 

```sh
python src/andes_gsea.py --emb embedding_file.csv --genelist embedding_gene_ids.txt --geneset gene_set_database.gmt --rankedlist ranked_genes.txt --out output_file.csv -n num_processor
```

