import re
import argparse


def adapt_output(fh_in, fh_adapt):
    for line in fh_in:
        if line.startswith('#'):
            continue
        else:
            wo_gaps = re.sub(' +', ' ', line)
            w_tabs = wo_gaps.replace(' ', '\t')
            w_tabs = re.split(r'\t+', w_tabs)
            important = w_tabs[0:1] + w_tabs[2:9] + w_tabs[11:15] + w_tabs[19:21]
            important = '\t'.join(important)
            fh_adapt.write('{}\n'.format(important))
    fh_adapt.close()
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fh_in', required=True, help='Fasta file')
    parser.add_argument('--fh_out', required=True, help='df output')
    args = parser.parse_args()
    tab_fh = open(args.fh_in)
    out_fh = open(args.fh_out, 'w')
    adapt_output(tab_fh, out_fh)
    
# fh = open(snakemake.input[0], 'r')
# fh_out = open(snakemake.output[0], 'w')
# adapt_output(fh, fh_out)
