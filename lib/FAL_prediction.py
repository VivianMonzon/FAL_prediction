import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import os


class FAL_prediction(object):

    def create_folder(self, analysisfolder, resultsfolder):
        # Create analysis folder:
        if os.path.exists('{}'.format(analysisfolder)):
            print('{} folder exists'.format(analysisfolder))
        else:
            os.system('mkdir {}'.format(analysisfolder))
            print('{} folder created'.format(analysisfolder))
        # Create results folder:
        if os.path.exists('{}'.format(resultsfolder)):
            print('{} folder exists'.format(resultsfolder))
        else:
            os.system('mkdir {}'.format(resultsfolder))
            print('{} folder created'.format(resultsfolder))

    def run_analysis(self, analysisfolder, resultsfolder, fasta_seqs, jobname,
                     treks, lipop, signalp, tmhmm):
        if jobname is None:
            jobname = fasta_seqs.split('.')[0]
        print('Collect feature:')
        os.system('lib/analyse_testing_data.sh {} {} {} {} {} {} {} {}'.format(
            fasta_seqs, analysisfolder, resultsfolder, jobname, treks, lipop,
            signalp, tmhmm))

    def run_random_forest(self,
                          analysisfolder,
                          resultsfolder, jobname):
        # train model:
        proteins = pd.read_csv('database/training_data.csv')
        proteins = proteins.drop_duplicates()
        proteins = proteins.drop('ID', axis=1)
        proteins = proteins.dropna()
        X_proteins = proteins.drop('Type', axis=1)
        y_proteins = proteins['Type']
        model = RandomForestClassifier(n_estimators=50, random_state=2,
                                       max_features=3)
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
        prediction_test_set('{}/{}_feature.csv'.format(
            analysisfolder, jobname), '{}/results_{}'.format(
                resultsfolder, jobname))
