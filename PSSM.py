# Import necessary libraries
from Bio import SeqIO, AlignIO, motifs
from Bio.Align.Applications import ClustalwCommandline
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Input Multi-FASTA File
fasta_file = "Phosphorylation.txt"

# Step 2: Perform Sequence Alignment using ClustalW
clustalw_exe = "/usr/bin/clustalw"  # Path to ClustalW executable
cline = ClustalwCommandline(clustalw_exe, infile=fasta_file)
stdout, stderr = cline()  # Run ClustalW

# The aligned file will be saved as "protein_sequences.aln" by ClustalW
aligned_file = "Phosphorylation.aln"
print("Alignment complete!")

# Step 3: Load the Aligned Sequences
alignment = AlignIO.read(aligned_file, "clustal")

# Step 4: Create a Motif and Generate a PSSM
aligned_sequences = [record.seq for record in alignment]
motif = motifs.create(aligned_sequences, alphabet='RHKDESTNQCUGPAILMFWYV-')

# View the frequency matrix (Position Frequency Matrix, PFM)
print("\nFrequency Matrix:")
print(motif.counts)
pd.DataFrame(motif.counts).to_csv('frequency_matrix_Phosphorylation.csv', index=False)

# Generate Position-Specific Scoring Matrix (PSSM)
pssm = motif.pssm
print("\nPosition-Specific Scoring Matrix (PSSM):")
for position, scores in pssm.items():
    print(f"Position {position}: {scores}")

# Step 5: Visualize the PSSM as a Heatmap
# Convert PSSM to a DataFrame for visualization
pssm_df = pd.DataFrame(pssm)
v = pssm_df.values

pssm_df.to_csv('pssm_matrix_Phosphorylation.csv', index=False)

# Plot the heatmap
plt.figure(figsize=(12, 8))
minval = np.nanmin(v[v != -np.inf])
maxval = np.nanmax(v[v != -np.inf])
sns.heatmap(pssm_df, cmap="viridis", annot=True, fmt=".2f", cbar=True, vmin=minval, vmax=maxval)
plt.title("Position-Specific Scoring Matrix (PSSM)")
plt.xlabel("Amino Acids")
plt.ylabel("Positions")
plt.savefig('Phosphorylation_HM.png', dpi=300)

#===============Acetylation, Sumoylation, ThiolOxidation, Methylation==============
#
# # Import necessary libraries
# from Bio import SeqIO, AlignIO, motifs
# from Bio.Align.Applications import ClustalwCommandline
# import seaborn as sns
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
#
# # Step 1: Input Multi-FASTA File
# fasta_file = "Methylation.txt"
#
# # Step 2: Perform Sequence Alignment using ClustalW
# clustalw_exe = "/usr/bin/clustalw"  # Path to ClustalW executable
# cline = ClustalwCommandline(clustalw_exe, infile=fasta_file)
# stdout, stderr = cline()  # Run ClustalW
#
# # The aligned file will be saved as "protein_sequences.aln" by ClustalW
# aligned_file = "Methylation.aln"
# print("Alignment complete!")
#
# # Step 3: Load the Aligned Sequences
# alignment = AlignIO.read(aligned_file, "clustal")
#
# # Step 4: Create a Motif and Generate a PSSM
# aligned_sequences = [record.seq for record in alignment]
# motif = motifs.create(aligned_sequences, alphabet='RHKDESTNQCUGPAILMFWYV-')
#
# # View the frequency matrix (Position Frequency Matrix, PFM)
# print("\nFrequency Matrix:")
# print(motif.counts)
# pd.DataFrame(motif.counts).to_csv('frequency_matrix_Methylation.csv', index=False)
#
# # Generate Position-Specific Scoring Matrix (PSSM)
# pssm = motif.pssm
# print("\nPosition-Specific Scoring Matrix (PSSM):")
# for position, scores in pssm.items():
#     print(f"Position {position}: {scores}")
#
# # Step 5: Visualize the PSSM as a Heatmap
# # Convert PSSM to a DataFrame for visualization
# pssm_df = pd.DataFrame(pssm)
# v = pssm_df.values
#
# pssm_df.to_csv('pssm_matrix_Methylation.csv', index=False)
#
# # Plot the heatmap
# plt.figure(figsize=(12, 8))
# minval = np.nanmin(v[v != -np.inf])
# maxval = np.nanmax(v[v != -np.inf])
# sns.heatmap(pssm_df, cmap="viridis", annot=True, fmt=".2f", cbar=True, vmin=minval, vmax=maxval)
# plt.title("Position-Specific Scoring Matrix (PSSM)")
# plt.xlabel("Amino Acids")
# plt.ylabel("Positions")
# plt.savefig('MethylationHM.png', dpi=300)
#



