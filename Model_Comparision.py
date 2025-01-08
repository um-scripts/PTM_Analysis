import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import 
matthews_corrcoef
from imblearn.metrics import specificity_score, sensitivity_score
import matplotlib.pyplot as plt


df1 = pd.read_csv(r"78UbiquitinationAdditional_SiteOutput.csv")
df2 = pd.read_csv(r"Musite_nofilter/MusiteDeep_Ubiquitination.csv")
df3 = pd.read_csv(r"210Merged.csv")

df2['UNIPROT'] = df2['ID'].apply(lambda x: x.split('|')[1])


o = []
w = []
s = []
mask1 = []

names = []
for i, m in df3.iterrows():
    seq = m['Sequence']
    flag_a = np.zeros(len(seq))
    flag_p = np.zeros(len(seq))
    b1 = df1[df1['UNIPROT']==m['UniprotID']]
    s1= b1['SITE'].dropna().apply(lambda x:int(x[1:])-1).to_numpy()
    b2 = df2[df2['UNIPROT']==m['UniprotID']]

    s2= b2['SITE'].dropna().apply(lambda x:int(x[1:])-1).to_numpy()

    if len(s1) > 0:
        flag_a[s1]=1
    if len(s2) > 0:
        flag_p[s2]=1
    o.append(flag_a)
    w.append(flag_p)
    s.append(seq)

    msk1 = np.zeros(len(seq))
    msk_idx = np.pad(flag_a.copy(), 5)
    for shift in range(-5, 6):
        msk1 += np.roll(msk_idx, shift)[5:-5]
    mask1.append(msk1 == 0)

curated_mask = np.concatenate(o)
pred_mask = np.concatenate(w)
SEQUENCE = np.array(list(''.join(s)))

print(SEQUENCE)
mask1 = np.concatenate(mask1)
mask2 = (SEQUENCE == 'K')
zero_mask = mask1 & mask2
one_mask = curated_mask == 1

assert np.sum(zero_mask & one_mask) == 0
# ======================================= Confusion Matrix ============================================

cf_matrix= confusion_matrix(curated_mask[zero_mask | one_mask], pred_mask[zero_mask | one_mask])
cf_mat= pd.DataFrame(cf_matrix, index = ['Curated_negative', 'Curated_positive'], columns=['Musite_negative', 'Musite_positive'])
accuracy= accuracy_score(curated_mask[zero_mask | one_mask], pred_mask[zero_mask | one_mask])
MCC= matthews_corrcoef(curated_mask[zero_mask | one_mask], pred_mask[zero_mask | one_mask])
sen_score = sensitivity_score(curated_mask[zero_mask | one_mask], pred_mask[zero_mask | one_mask], labels= ['class 0', 'class 1'])
spec_score = specificity_score(curated_mask[zero_mask | one_mask], pred_mask[zero_mask | one_mask], labels= ['class 0','class 1'])
sns.heatmap(cf_mat, annot=True, fmt=".0f", cmap=ListedColormap(['white']), linewidths=1,linecolor='black', cbar=False)
plt.title("MusiteDeep")
plt.savefig('Musite/confusionmatrix_Musite_Ubiquitination_RESULT.png')
plt.close()

print(cf_matrix)
print(accuracy, MCC, sen_score, spec_score)

#========================bootstrap=============================

actual_index = np.arange(len(SEQUENCE))
all_accuracy= []
all_MCC= []
all_sensitivity= []
all_specificity= []

sample_size = int(0.9 * min(np.sum(one_mask), np.sum(zero_mask)))
for n in range(1000):
    ya_positive = np.random.choice(actual_index[one_mask], size= sample_size)
    ya_negetive = np.random.choice(actual_index[zero_mask], size= sample_size)
    curated = np.concatenate([curated_mask[ya_positive], curated_mask[ya_negetive]])
    predicted = np.concatenate([pred_mask[ya_positive], pred_mask[ya_negetive]])
    acc_all = accuracy_score(curated, predicted)
    all_accuracy.append(acc_all)
    MCC_all = matthews_corrcoef(curated, predicted)
    all_MCC.append(MCC_all)
    sensitivity_all = sensitivity_score(curated, predicted, labels= ['class 0','class 1'])
    all_sensitivity.append(sensitivity_all)
    specificity_all = specificity_score(curated, predicted, labels= ['class 0','class 1'])
    all_specificity.append(specificity_all)
    print('Accuracy:', acc_all, 'sensitivity:', sensitivity_all, 'Specificity:' ,specificity_all)
median_accuracy= np.median(all_accuracy)
median_specificity= np.median(all_specificity)
median_sensitivity= np.median(all_sensitivity)
median_MCC= np.median(all_MCC)
colors= ['peachpuff', 'cyan', 'lightgreen', 'yellow']
bplot= plt.boxplot([all_MCC, all_accuracy, all_sensitivity, all_specificity], labels=['MCC','Accuracy','Sensitivity','Specificity'],patch_artist=True)

for patch,color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

plt.ylim([0,1.1])


print(median_accuracy)
print(median_MCC)
print(median_specificity)
print(median_sensitivity)
plt.savefig('Musite/Musite_Ubiquitination_result.png')
plt.close()

