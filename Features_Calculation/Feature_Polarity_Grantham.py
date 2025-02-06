import pandas as pd
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import os

grantham_polarity = {
    'A': 8.1, 'R': 10.5, 'N': 11.6, 'D': 13.0, 'C': 5.5, 
    'Q': 10.5, 'E': 12.3, 'G': 9.0, 'H': 10.4, 'I': 5.2, 
    'L': 4.9, 'K': 11.3, 'M': 5.7, 'F': 5.2, 'P': 8.0, 
    'S': 9.2, 'T': 8.6, 'W': 5.4, 'Y': 6.2, 'V': 5.9
}

def calculate_polarity_scores(sequence, window_size=7):
    if window_size > len(sequence):
        window_size = len(sequence)
    
    polarity_values = [grantham_polarity.get(aa.upper(), 0) for aa in sequence]
    
    scores = []
    half_window = window_size // 2
    
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        window = polarity_values[start:end]
        scores.append(sum(window) / len(window))
    
    return scores

def analyze_fasta(input_file, output_dir='polarity_results'):
    os.makedirs(output_dir, exist_ok=True)
    
    all_sequence_data = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        polarity_scores = calculate_polarity_scores(str(record.seq))
        row_data = [record.id] + polarity_scores
        all_sequence_data.append(row_data)
        plt.figure(figsize=(12, 4))
        plt.plot(polarity_scores)
        avg_polarity = sum(polarity_scores) / len(polarity_scores)
        plt.axhline(y=avg_polarity, color='r', linestyle='-', label=f'Average Polarity ({avg_polarity:.2f})')
        plt.axhline(y=10.0, color='g', linestyle='--', label='High Polarity Threshold')
        plt.axhline(y=6.0, color='b', linestyle='--', label='Low Polarity Threshold')
        plt.title(f"Polarity Profile (Grantham) for Sequence: {record.id}")
        plt.xlabel("Position")
        plt.ylabel("Polarity Score")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plot_filename = os.path.join(output_dir, f"{record.id}_polarity_plot.png")
        plt.savefig(plot_filename)
        plt.close()
    
    column_names = ['Sequence_ID'] + [f'Pos_{i+1}_Score' for i in range(len(polarity_scores))]
    df = pd.DataFrame(all_sequence_data, columns=column_names)
    csv_filename = os.path.join(output_dir, 'polarity_scores.csv')
    df.to_csv(csv_filename, index=False)
    
    print(f"Analysis complete. Results in {output_dir}:")
    print(f"- CSV file: {csv_filename}")
    print(f"- Total sequences: {len(all_sequence_data)}")
    print(f"- Plots saved in: {output_dir}")

def main():
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python script.py input.fasta [output_directory]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'polarity_results'
    
    try:
        analyze_fasta(input_file, output_dir)
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()
