from Bio import SeqIO
import matplotlib.pyplot as plt
import csv
import numpy as np

# Kyte-Doolittle hydrophobicity scale
kd_scale = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def calculate_hydrophobicity(sequence, window_size=7):
    fixed_length = 41
    
    if len(sequence) > fixed_length:
        sequence = sequence[:fixed_length]
    else:
        sequence = sequence.ljust(fixed_length, '-')
    
    if window_size > len(sequence):
        window_size = len(sequence)

    hydro_values = [kd_scale.get(aa, 0) for aa in sequence]
    
    scores = []
    half_window = window_size // 2
    
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        window = hydro_values[start:end]
        scores.append(sum(window) / len(window))
    
    return scores

def analyze_fasta_file(input_file, window_size=7):
    all_hydro_scores = []
    all_seq_ids = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        seq_id = record.id
        
        hydrophobicity_scores = calculate_hydrophobicity(str(record.seq), window_size)
        
        all_hydro_scores.append(hydrophobicity_scores)
        all_seq_ids.append(seq_id)
        
        plt.figure(figsize=(12, 4))
        plt.plot(hydrophobicity_scores)
        plt.axhline(y=0, color='r', linestyle='-')
        plt.title(f"Hydrophobicity Profile for Sequence: {seq_id}")
        plt.xlabel("Position")
        plt.ylabel("Hydrophobicity Score")
        plt.savefig(f"{seq_id}_hydrophobicity_plot.png")
        plt.close()
    
    csv_filename = "hydrophobicity_scores.csv"
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        header = ['Sequence_ID'] + [f'Pos_{i+1}_Score' for i in range(41)]
        writer.writerow(header)
        
        for seq_id, scores in zip(all_seq_ids, all_hydro_scores):
            row = [seq_id] + [f'{score:.3f}' for score in scores]
            writer.writerow(row)

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python script.py input.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    try:
        analyze_fasta_file(input_file)
        print("Analysis complete. Plots and CSV file have been saved.")
    except Exception as e:
        print(f"Error: {str(e)}")
