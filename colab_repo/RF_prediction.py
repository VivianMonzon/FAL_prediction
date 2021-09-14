import pandas as pd
from sklearn.ensemble import RandomForestClassifier


class ML_prediction():
     
    def run_random_forest(features, training_set,
                          resultsfolder, output_fh):
        proteins = pd.read_csv(training_set)
        proteins = proteins.drop_duplicates()
        proteins = proteins.drop('ID', axis=1)
        proteins = proteins.dropna()
        X_proteins = proteins.drop('Type', axis=1)
        y_proteins = proteins['Type']
        model = RandomForestClassifier(n_estimators=50,
                                       random_state=2, max_features=3)
        model.fit(X_proteins, y_proteins)
        test_set = pd.read_csv(features)
        test_set = test_set.dropna()
        test_set_wo_id = test_set.drop('ID', axis=1)
        y_model = model.predict(test_set_wo_id)
        probability_1 = model.predict_proba(test_set_wo_id)[:, 1]
        test_set['prediction'] = y_model
        test_set['probability_Adh'] = probability_1
        test_set = test_set[['ID', 'prediction', 'probability_Adh',
                             'length', 'Hydro_portion', 'Charge_portion', 
                             'rel_entropy', 'R', 'H', 'K', 'D', 'E', 'S', 'T', 
                             'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 
                             'F', 'Y', 'W']]
        test_set.to_csv('{}/{}.tsv'.format(resultsfolder,
                                           output_fh), index=False, sep='\t')
