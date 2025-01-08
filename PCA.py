import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import pandas as pd
import numpy as np

DATA = pd.read_csv('features220.csv')
print(DATA)
n_comp = 2
pca = PCA(n_components=n_comp)

F = pca.fit_transform(DATA)
for i in range(n_comp):
	feature_weights = np.abs(pca.components_[i])
	feature_importance = feature_weights * 100 / np.sum(feature_weights)
	idx = np.argsort(feature_importance)[-15:]
	xl = pca.feature_names_in_[idx]
	yl = feature_importance[idx]

	plt.title(f"Component {i}: Explained Variance: {pca.explained_variance_ratio_[i]*100:.3f}%")
	plt.barh(xl, yl, edgecolor='k')
	plt.gca().xaxis.grid()
	plt.xlabel('Feature Importance (%)')
	plt.tight_layout()
	plt.show()

#kmeans = KMeans(n_clusters=3, n_init='auto')
#cl = kmeans.fit_predict(F)

#N = pd.read_csv('220_identifiers.csv')
#for i in range(n_comp):
#    N[f'COMP_{i}'] = F[:, i]
#N['CLUSTER'] = cl
#N.to_csv('clusters.csv', index=False)


#------------------PCA PLOT-----------------
#import pandas as pd
#import matplotlib.pyplot as plt

#T = pd.read_csv('clusters.csv')
#for cl in range(3):
#    d = T[T['CLUSTER'] == cl]
#    plt.scatter(d['COMP_0'], d['COMP_1'], edgecolors='k', label=f"cluster {cl}");
#plt.legend()
#plt.xlabel('Principal Component 0')
#plt.ylabel('Principal Component 1')
#plt.savefig('fig.png', dpi=500)
