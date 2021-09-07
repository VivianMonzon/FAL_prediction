import argparse
from lib.controller import Controller


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="commands")

    # Arguments for feature run analysis:
    run_analysis = subparsers.add_parser("collect_feature", 
                                         help='run analysis for fasta seqs to get feature')
    run_analysis.add_argument('--analysisfolder', default='analysis',
                              help='give path to analysis folder, ' +
                              'if not set analysis folder will be created')
    run_analysis.add_argument('--resultsfolder', default='results',
                              help='give path to results folder, ' +
                              'if not set results folder will be created')
    run_analysis.add_argument('--fasta_seqs', required=True, 
                              help='path to input fasta sequences')
    run_analysis.set_defaults(func=collect_feature)

    # Arguments for RandomForest prediction:
    ML_predict = subparsers.add_parser("ML_prediction", 
                                       help='ML prediction for FAL proteins')
    ML_predict.add_argument('--analysisfolder', default='analysis',
                            help='give path to analysis folder, ' +
                            'if not set analysis folder will be created')
    ML_predict.add_argument('--resultsfolder', default='results',
                            help='give path to results folder, ' +
                            'if not set results folder will be created')
    ML_predict.set_defaults(func=ML_prediction)

    args = parser.parse_args()
    if "func" in dir(args):
        args.func(args)
    else:
        parser.print_help()
	

def collect_feature(args):
    run_feature_analysis = Controller(args)
    run_feature_analysis.collect_feature()


def ML_prediction(args):
    run_ML_prediction = Controller(args)
    run_ML_prediction.ML_prediction()


main()
