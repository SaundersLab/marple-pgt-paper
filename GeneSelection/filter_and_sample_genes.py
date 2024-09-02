from read_gff import read_gff
import pandas as pd
from Bio import SeqIO
import glob
from pathlib import Path
import sys
from collections import defaultdict
import numpy as np
import matplotlib
# use agg backend because tkinter isn't on HPC
matplotlib.use('agg')
import matplotlib.pyplot as plt

# Load GFF data about CDS
def load_cds_data(gff_path):    
    gff = read_gff(gff_path)
    gff = gff.reset_index()
    genes = gff[gff['type'] == 'gene'].set_index('ID')
    cds_groups = gff[gff['type'] == 'CDS'].groupby('ID') # changed groupby to set_index
    cds = pd.DataFrame({
        'index': cds_groups.index.first(),
        'seqid': cds_groups.seqid.first(),
        'product': cds_groups.product.first(),
        'protein_id': cds_groups.protein_id.first(),
        'cds_id': cds_groups.ID.first(),
        'gene_id': [[a for a in l if a[:4] == 'gene'][0] for l in cds_groups.Ancestors.first()],
        'cds_length': cds_groups.length.sum(),
    })
    cds['gene_length'] = [genes.loc[ID].length for ID in cds.gene_id]
    cds['gene_name'] = [genes.loc[ID].Name for ID in cds.gene_id]
    cds = cds.reset_index().drop(columns='ID')
    cds = cds.sort_values(by='index')
    # Create a name for each CDS variant of each gene, based on the order it appeared in the gene, e.g. 
    # the 3 CDS variants for the gene PGTG_12162 will get the names PGTG_12162_1, PGTG_12162_2, and PGTG_12162_3.
    # For these names to match up with the consensus files, we must use the same GFF file
    gene_variant_names = []
    for gene_name, prev_gene_name in zip(cds.gene_name, ('',) + tuple(cds.gene_name[:-1])):
        gene_count = 1 if gene_name != prev_gene_name else gene_count + 1
        gene_variant_names.append(gene_name + '_' + str(gene_count))
    cds['gene_variant_name'] = gene_variant_names
    # Keep track of how many variants each gene has
    gene_id_to_n_variants = cds.groupby('gene_id').size().to_dict() # changed groupby to set_index
    cds['n_variants']  = [gene_id_to_n_variants[gene_id] for gene_id in cds.gene_id]
    return cds

print(sys.argv)
gff_path = sys.argv[1]
alignments_out_dir = sys.argv[2]
reports_out_dir = sys.argv[3]
consensus_dir = sys.argv[4]
outgroup_consensus_dir = sys.argv[5]

if not Path(reports_out_dir).exists():
    Path(reports_out_dir).mkdir(parents=True)
if not Path(alignments_out_dir).exists():
    Path(alignments_out_dir).mkdir(parents=True)

fasta_paths = glob.glob(f'{consensus_dir}/*.fa')
outgroup_fasta_paths = glob.glob(f'{outgroup_consensus_dir}/*.fa')

isolates = [path.split('/')[-1].rsplit('_', 2)[0] for path in fasta_paths]
min_p_cov = .8

cds = load_cds_data(gff_path)

# Filter gene and CDS length
cds = cds[(cds.cds_length >= 1000) & (cds.gene_length <= 2000)]

# Calculate coverage
n_isolates = len(fasta_paths)

def p_cov(seq):
    return 1 - seq.count('?') / len(seq)

# For each isolate, create a column with the coverage for each gene (1 row per gene)
for path, isolate in zip(fasta_paths, isolates):
    fasta_dict = SeqIO.index(path, 'fasta')
    cds[isolate] = [p_cov(fasta_dict[gv].seq) for gv in cds.gene_variant_name]
    
# For each gene, find the percentage of isolates that have at least 
# min_p_cov proportion of the gene's sites covered
site_covs = cds[cds.columns[-n_isolates:]].values
p_isolates_covered = [(covs > min_p_cov).sum() / n_isolates for covs in site_covs]
cds.insert(11, 'p_isolates_covered', p_isolates_covered)

# Filter genes based on percent coverage across all isolates
cds = cds[cds.p_isolates_covered >= min_p_cov]

# Calculate SNPs per base - only use every 3rd base in the calculation
site_stats = {
    gv: {
        # which bases (excluding '?') are present at this site across isolates
        'unique_bases': [set() for _ in range(2, l, 3)],
        # proportion of isolates with a non-missing (i.e. not '?') base at this site
        'p_cov': [0 for _ in range(2, l, 3)]
    }
    for gv, l in zip(cds.gene_variant_name, cds.cds_length)
}

# For every 3rd base in every gene, find the different bases at that site across
# all isolates, as well as the proportion of isolates which are covered at that site
for path, isolate in zip(fasta_paths, isolates):
    fasta_dict = SeqIO.index(path, 'fasta')
    for gv in cds.gene_variant_name:
        seq = fasta_dict[gv].seq
        for i, base in enumerate(seq[2::3]):
            if base != '?':
                site_stats[gv]['unique_bases'][i].add(base)  
                site_stats[gv]['p_cov'][i] += 1 / n_isolates
                
snps_per_base = {gv: 0 for gv in site_stats}
for gv, stats in site_stats.items():
    for unique_bases, p_cov in zip(stats['unique_bases'], stats['p_cov']):
        if len(unique_bases) > 1 and p_cov >= min_p_cov:
            snps_per_base[gv] += 1 / len(stats['unique_bases'])
            
cds.insert(11, 'snps_per_base', [snps_per_base[gv] for gv in cds.gene_variant_name])

# Estimate the probability of being included in a sample
# Longer genes are less likely to be chosen, and the average probability 
# of a gene being chosen is just under 60%
cds.insert(12, 'p_chosen', .6 - .2 * (cds.gene_length - 1500) / 1000)

# Create random samples of the genes that have qualifying SNPs per base for a range of thresholds
# e.g. "0.1_4" is a column of flags which are True if the gene has at least 0.1 SNPs per base
# and was randomly selected to be in the 4th sample at that threshold
thresholds = np.arange(.09, .19, .005)
n_samples = 6
n = cds.shape[0]
np.random.seed(0)
for thresh in thresholds:
    qualifies = cds.snps_per_base >= thresh
    cds[f'{round(thresh, 3)}'] = qualifies
    for i in range(n_samples):
        is_randomly_chosen = cds.p_chosen > np.random.uniform(0, 1, n)
        cds[f'{round(thresh, 3)}_{i + 1}'] = qualifies & is_randomly_chosen
        
cds.to_csv(f'{reports_out_dir}/sampled_cds.csv', index=False)

# Create alignments from the samples
for col in cds.columns[-(n_samples * (len(thresholds) + 1)):]:
    gene_variant_names = cds[cds[col]].gene_variant_name.values
    alignment = ''
    for path in (fasta_paths + outgroup_fasta_paths):
        isolate = path.split('/')[-1].split('_')[0]
        fasta_dict = SeqIO.index(path, 'fasta')
        alignment += f'>{isolate}\n'
        for gv in gene_variant_names:
            alignment += str(fasta_dict[gv].seq[2::3])
        alignment += '\n'
    with open(f'{alignments_out_dir}/{col}.fa', 'w') as f:
        f.write(alignment)
        
# Create plots to illustrate properties of the samples
n_qualifying_genes = []
mean_gene_lengths = []
mean_cds_lengths = []
mean_n_variants = []

fig, axes = plt.subplots(2, 2, figsize=(20, 20))
for threshold in thresholds:
    n_qualifying_genes.append(cds[cds.snps_per_base >= threshold].shape[0])
    mean_gene_lengths.append(cds[cds.snps_per_base >= threshold].gene_length.mean())
    mean_cds_lengths.append(cds[cds.snps_per_base >= threshold].cds_length.mean())
    mean_n_variants.append(cds[cds.snps_per_base >= threshold].n_variants.mean())

ax = axes[0][0]
ax.plot(thresholds, n_qualifying_genes, marker='o', color='k', label='Qualifying Genes', markersize=10)
for i in range(1, n_samples + 1):
    ax.scatter(thresholds, [cds[f'{round(thresh, 3)}_{i}'].sum() for thresh in thresholds], label=f'Sample {i}', marker='o', alpha=0.5, s=100)
ax.grid(); ax.legend(); ax.set_xlabel('SNPs per Base')
ax.set_ylabel('Genes')

ax = axes[0][1]
ax.plot(thresholds, mean_gene_lengths, marker='o', color='k', label='Qualifying Genes', markersize=10)
for i in range(1, n_samples + 1):
    ax.scatter(thresholds, [cds[cds[f'{round(thresh, 3)}_{i}']].gene_length.mean() for thresh in thresholds], label=f'Sample {i}', marker='o', alpha=0.5, s=100)
ax.grid(); ax.legend(); ax.set_xlabel('SNPs per Base')
ax.set_ylabel('Mean Gene Length')

ax = axes[1][0]
ax.plot(thresholds, mean_cds_lengths, marker='o', color='k', label='Qualifying Genes', markersize=10)
for i in range(1, n_samples + 1):
    ax.scatter(thresholds, [cds[cds[f'{round(thresh, 3)}_{i}']].cds_length.mean() for thresh in thresholds], label=f'Sample {i}', marker='o', alpha=0.5, s=100)
ax.grid(); ax.legend(); ax.set_xlabel('SNPs per Base')
ax.set_ylabel('Mean CDS Length')

ax = axes[1][1]
ax.plot(thresholds, mean_n_variants, marker='o', color='k', label='Qualifying Genes', markersize=10)
for i in range(1, n_samples + 1):
    ax.scatter(thresholds, [cds[cds[f'{round(thresh, 3)}_{i}']].n_variants.mean() for thresh in thresholds], label=f'Sample {i}', marker='o', alpha=0.5, s=100)
ax.grid(); ax.legend(); ax.set_xlabel('SNPs per Base')
ax.set_ylabel('Mean Number of CDS Variants')

plt.savefig(f'{reports_out_dir}/sampled_cds_plots.svg')
