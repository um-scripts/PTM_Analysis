from Bio import SeqIO, AlignIO, motifs
from Bio.Align.Applications import ClustalwCommandline
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fasta_file = "Phosphorylation.txt"

#ClustalW
clustalw_exe = "/usr/bin/clustalw"  # Path to ClustalW executable
cline = ClustalwCommandline(clustalw_exe, infile=fasta_file)
stdout, stderr = cline()  # Run ClustalW
aligned_file = "Phosphorylation.aln"
alignment = AlignIO.read(aligned_file, "clustal")

#PSSM
aligned_sequences = [record.seq for record in alignment]
motif = motifs.create(aligned_sequences, alphabet='RHKDESTNQCUGPAILMFWYV-')
print("\nFrequency Matrix:")
print(motif.counts)
pd.DataFrame(motif.counts).to_csv('frequency_matrix_Phosphorylation.csv', index=False)
pssm = motif.pssm
print("\nPosition-Specific Scoring Matrix (PSSM):")
for position, scores in pssm.items():
    print(f"Position {position}: {scores}")
pssm_df = pd.DataFrame(pssm)
v = pssm_df.values

pssm_df.to_csv('pssm_matrix_Phosphorylation.csv', index=False)

plt.figure(figsize=(12, 8))
minval = np.nanmin(v[v != -np.inf])
maxval = np.nanmax(v[v != -np.inf])
sns.heatmap(pssm_df, cmap="viridis", annot=True, fmt=".2f", cbar=True, vmin=minval, vmax=maxval)
plt.title("Position-Specific Scoring Matrix (PSSM)")
plt.xlabel("Amino Acids")
plt.ylabel("Positions")
plt.savefig('Phosphorylation_HM.png', dpi=300)


