import scanpy as sc

x = sc.read_10x_mtx('../data/mm10/')
# z = sc.read_10x_mtx('../data/mm10/matrix.mtx', var_names='../data/mm10/genes.tsv'[1])
print()