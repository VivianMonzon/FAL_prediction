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
    FA_prediction.add_argument('--lipop_dir', required=True,
                               help='path to folder containing LipoP tool.'
                               'Example: /Users/vmonzon/Downloads/software/LipoP1.0a')
    FA_prediction.add_argument('--signalp_dir', required=True,
                               help='path to folder containing SignalP tool. '
                               'Example: /Users/vmonzon/Downloads/software/signalp-5.0/bin')
    FA_prediction.add_argument('--tmhmm_dir', required=True,
                               help='path to folder containing TMHMM tool. '
                               'Example: /Users/vmonzon/Downloads/software/tmhmm-2.0c/bin')
    FA_prediction.add_argument('--iupred_dir', required=True,
                               help='path to folder containing iupred2a.py file'
                               'and corresponding data folder. Example: '
                               '/Users/vmonzon/Downloads/software')
    FA_prediction.set_defaults(func=predict)
    
    args = parser.parse_args()
    if "func" in dir(args):
        args.func(args)
    else:
        parser.print_help()


def predict(args):
    FA_prediction_analysis = Controller(args)
    FA_prediction_analysis.predict()


main()
