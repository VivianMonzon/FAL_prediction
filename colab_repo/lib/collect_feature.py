import os


class Collect_feature(object):

    def run_analysis(self, analysisfolder, resultsfolder, fasta_seqs):
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
        # Collect feature:
        os.system('lib/analyse_testing_data.sh {} {} {}'.format(fasta_seqs,
                                                                analysisfolder, resultsfolder))
