import pandas as pd
from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import accuracy_score
# from sklearn.tree import export_graphviz
# import matplotlib.pyplot as plt
# import seaborn as sns; sns.set(font_scale=1.2)


def prediction_test_set(fh_in, fh_out):
    test_set = pd.read_csv(fh_in)
    test_set = test_set.dropna()
    test_set_wo_id = test_set.drop('ID', axis=1)
    y_model = model.predict(test_set_wo_id)
    probability_1 = model.predict_proba(test_set_wo_id)[:,1]
    test_set['prediction'] = y_model
    test_set['probability_Adh'] = probability_1
    test_set.to_csv(fh_out, index=False, sep='\t')


proteins = pd.read_csv(snakemake.input[0])
proteins = proteins.drop_duplicates()
proteins = proteins.drop('ID', axis=1)
proteins = proteins.dropna()

X_proteins = proteins.drop('Type', axis=1)
y_proteins = proteins['Type']

model = RandomForestClassifier(n_estimators=50, random_state=2, max_features=3)
model.fit(X_proteins, y_proteins)

prediction_test_set(snakemake.input[1], snakemake.output[0])
