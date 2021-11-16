import argparse
from lib.controller import Controller


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="commands")

    FA_prediction = subparsers.add_parser("predict",
                                          help='FA-like prediction analysis')
    FA_prediction.add_argument('--analysisfolder', default='analysis',
                               help='give path to analysis folder, ' +
                               'if not set analysis folder will be created')
    FA_prediction.add_argument('--resultsfolder', default='results',
                               help='give path to results folder, ' +
                               'if not set results folder will be created')
    FA_prediction.add_argument('--fasta_seqs', required=True,
                               help='path to input fasta sequences')
    FA_prediction.add_argument('--jobname', required=False,
                               help='If none given uses fasta file name')
    FA_prediction.add_argument('--treks_dir', required=True,
                               help='path to folder containing T-REKS tool. ' 
                               'Example: /Users/vmonzon/Downloads')
    FA_prediction.add_argument('--lipop_dir', required=False,
                               help='path to folder containing LipoP tool.'
                               'If none given, excepts to find it in same '
                               'directory as T-REKS. Example: '
                               '/Users/vmonzon/Downloads/software/LipoP1.0a')
    FA_prediction.add_argument('--signalp_dir', required=False,
                               help='path to folder containing SignalP tool. '
                               'If none given, excepts to find it in same '
                               'directory as T-REKS. Example: '
                               '/Users/vmonzon/Downloads/software/signalp-5.0/bin')
    FA_prediction.add_argument('--tmhmm_dir', required=False,
                               help='path to folder containing TMHMM tool. '
                               'If none given, excepts to find it in same '
                               'directory as T-REKS. Example: '
                               '/Users/vmonzon/Downloads/software/tmhmm-2.0c/bin')
    FA_prediction.set_defaults(func=predict)
    
    args = parser.parse_args()
    if "func" in dir(args):
        args.func(args)
    else:
        parser.print_help()


def predict(args):
    FA_prediction_analysis = Controller(args)
    FA_prediction_analysis.predict()

    
# def create(args):
#     create_folder_controller = Controller(args)
#     create_folder_controller.create()


# def collect_feature(args):
#     run_feature_analysis = Controller(args)
#     run_feature_analysis.collect_feature()


# def ML_prediction(args):
#     run_ML_prediction = Controller(args)
#     run_ML_prediction.ML_prediction()


main()
