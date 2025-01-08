from sklearn.preprocessing import OneHotEncoder
import pandas as pd
import numpy as np


D = pd.read_csv(r"C:\Users\dassa\Downloads\Glycosylation3Di.csv")

DF = D[D['3Di'].apply(lambda s: len(s) == 41)]
F = np.vstack(DF['3Di'].apply(lambda seq: list(seq)).to_numpy())
ohencoder = OneHotEncoder(categories=[list('ACDEFGHIKLMNPQRSTVWY')] * 41, sparse_output=False)
XX = ohencoder.fit_transform(F)
print(XX.shape)
