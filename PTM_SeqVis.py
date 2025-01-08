import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

ptm_colors = {
    'Phosphorylation': ('#FF0000', 'o'),     # Red
    'Methylation': ('#0000FF', 's'),         # Blue
    'Acetylation': ('#00FF00', '^'),         # Green
    'Sumoylation': ('#FFA500', '<'),         # Orange
    'Ubiquitination': ('#800080', '>'),      # Purple
    'Palmitoylation': ('#FF8DA1', 'v'),      # Gold
    'Glycosylation': ('#00FFFF', 'D'),        # Cyan
    'Monomethylation': ('#BFFF00', 'p')
}
        
SEQ = pd.read_csv('48_NRHuman.csv')
df = pd.read_csv("NR1F1.csv")
PROT = df.groupby('UniProt_ID')
for uid in PROT.groups:
    print(uid)
    sequence = SEQ[SEQ['UniProt_ID'] == uid]['Sequence'].values[0]
    NPOS_PERLINE = 45
    SPACING = 2.2
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111)
    ax.set_xlim([-1, 130])
    ax.set_ylim([-1, 28])
    ax.invert_yaxis()
    
    print(sequence)
    for i, aa in enumerate(sequence):
        xpos = (i % NPOS_PERLINE) * SPACING
        ypos = (i // NPOS_PERLINE) * SPACING
        ax.text(xpos, ypos, aa, ha='center', va='center', weight="bold", fontsize=12, zorder=2)
    
    
    DD = PROT.get_group(uid)
    legend_elements = []
    for ptm, (color, marker) in ptm_colors.items():
        mod_positions = DD[DD['PTM'] == ptm]['SITE'].values
        # Convert to Zero-based Indexing
        mod_positions = np.unique(np.array([pos - 1 for pos in mod_positions]).astype(int))
        print(ptm, mod_positions)
        
        diagram_positions_x, diagram_positions_y = [], []
        for pos in mod_positions:
            xpos = (pos % NPOS_PERLINE) * SPACING
            ypos = (pos // NPOS_PERLINE) * SPACING
            print(pos, xpos, ypos)
            diagram_positions_x.append(xpos)
            diagram_positions_y.append(ypos)
        
        lh = ax.scatter(diagram_positions_x, diagram_positions_y, color=color, s=250, alpha=0.4, zorder=1, label=ptm, marker=marker, edgecolors='k')
        legend_elements.append(lh)
    
    ax.set_title(f'Protein Modifications for {uid}')
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', title='Modifications', bbox_to_anchor=(0.45, 0.4, 0.5, 0.5), labelspacing=1.1, borderpad=1.1, edgecolor='k')
    plt.subplots_adjust(top=0.9, right=0.95, left=0.05, bottom=0.05)
    ax.axis('off')
    # plt.show()
    plt.savefig(f'{uid}_modifications.png', dpi=500, bbox_inches='tight')
    plt.close()

