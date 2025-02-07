# import requests
#
# with open('51_notPredicted.txt') as f:
#     ids = f.readline().strip('\n').split(',')
# for iid in set(ids):
#     print(iid)
#     result=requests.get(f'http://iupred2a.elte.hu/iupred2a/short/{iid}')
#     res=result.text
#     res=res.replace('<pre>','')
#     res=res.replace('</pre>','')
#     with open(f'51_notPredicted_short_results/{iid}.tsv','w') as file:
#         file.write(res)
#         file.close()

# ============================================extract IDR values for Curated predicted regions=======================
# import matplotlib.pyplot as plt
# import pandas as pd
# import seaborn as sns
#
# # EXTRACT IDRs for curated sites
# df = pd.read_csv('51_unmappedSites.csv')
# df['pos'] = df['SITE'].apply(lambda x: int(x[1:]))
#
# for i, r in df.iterrows():
#     dd = pd.read_csv(f'51_notPredicted_short_results/{r["UNIPROT"]}.tsv', comment='#', sep=r'\s+', header=None, names=['col1', 'col2', 'col3'])
#     df.loc[i, 'V1'] = dd.loc[dd['col1'] == r['pos'], 'col2'].iloc[0]
#     df.loc[i, 'V2'] = dd.loc[dd['col1'] == r['pos'], 'col3'].iloc[0]
#
# df.to_csv('IDRs_51notpredicted_phospho_short.csv', index=False)

# ================plot IDRs=======================
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd

cmap = cm.get_cmap('RdPu', 5)
norm = plt.Normalize(0, 1)

fig = plt.figure(figsize=(8, 10))
gs = GridSpec(1, 2, width_ratios=[20, 1], wspace=0.01)
ax = fig.add_subplot(gs[0])
lax = fig.add_subplot(gs[1])
# lax.axis('off')

df = pd.read_csv('IDRs_51notpredicted_phospho_short.csv')
# sns.barplot(df, x="CuratedList", y="IUPredScore")
# proteins = np.array(df['UNIPROT'].unique())
ax.scatter(df['pos'], df['UNIPROT'], c=cmap(norm(df['IUPredScore'])), edgecolors='k')
cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=lax, format=lambda x, pos: '{:.1f}'.format(x), orientation='vertical')
cbar.set_label('IUPredScore', rotation=90, fontsize=12)

plt.subplots_adjust(left=0.12, right=0.9, bottom=0.06, top=0.95)
plt.savefig('IDRs_51notpredicted_phospho_short.png')

