from lib.FAL_prediction import FAL_prediction


class Controller(object):

    def __init__(self, args):
        self._args = args

    def predict(self):
        FA_prediction = FAL_prediction()
        FA_prediction.create_folder(
            self._args.analysisfolder,
            self._args.resultsfolder)
        FA_prediction.run_analysis(
            self._args.analysisfolder,
            self._args.resultsfolder,
            self._args.fasta_seqs,
            self._args.jobname,
            self._args.treks_dir,
            self._args.lipop_dir,
            self._args.signalp_dir,
            self._args.tmhmm_dir,
            self._args.iupred_dir)
        FA_prediction.run_random_forest(
            self._args.analysisfolder,
            self._args.resultsfolder,
            self._args.jobname)
