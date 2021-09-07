import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import seaborn as sns
sns.set(font_scale=1.2)


class ML_prediction(object):
        
    def run_random_forest(self, 
                          analysisfolder, 
                          resultsfolder):
        # train model:
        proteins = pd.read_csv('../resources/training_data.csv')
        proteins = proteins.drop_duplicates()
        proteins = proteins.drop('ID', axis=1)
        proteins = proteins.dropna()
        X_proteins = proteins.drop('Type', axis=1)
        y_proteins = proteins['Type']
        model = RandomForestClassifier(n_estimators=50,
                                       random_state=2, max_features=3)
        model.fit(X_proteins, y_proteins)

        def prediction_test_set(fh_in, fh_out):
            test_set = pd.read_csv(fh_in)
            test_set = test_set.dropna()
            test_set_wo_id = test_set.drop('ID', axis=1)
            y_model = model.predict(test_set_wo_id)
            probability_1 = model.predict_proba(test_set_wo_id)[:, 1]
            test_set['prediction'] = y_model
            test_set['probability_Adh'] = probability_1
            test_set.to_csv('{}.tsv'.format(fh_out), index=False, sep='\t')
        
        # predict:
        prediction_test_set('{}/query_seq_feature.csv'.format(analysisfolder), 
                            '{}/results'.format(resultsfolder))
