from lib.collect_feature import Collect_feature
from lib.random_forest_prediction import ML_prediction


class Controller(object):

    def __init__(self, args):
        self._args = args

    def collect_feature(self):
        run_feature_analysis = Collect_feature()
        run_feature_analysis.run_analysis(
            self._args.analysisfolder, 
            self._args.resultsfolder, 
            self._args.fasta_seqs)

    def ML_prediction(self):
        run_ML_prediction = ML_prediction()
        run_ML_prediction.run_random_forest(
            self._args.analysisfolder,
            self._args.resultsfolder)
