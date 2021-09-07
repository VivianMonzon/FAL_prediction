import pandas as pd
import argparse
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument('--domtblout', required=True, help="hmmer domtblout")
parser.add_argument('--out_readable', required=True, help="readable output")
parser.add_argument('--out_important', required=True, help="only important columns")

args = parser.parse_args()

fh = open(args.domtblout, 'r')
wo_tab_file = open(args.out_readable, 'w')
for line in fh:
    if line.startswith('#'):
        continue
    else:
        wo_gaps = re.sub(' +', ' ', line)
        w_tabs = wo_gaps.replace(' ', '\t')
        wo_tab_file.write(w_tabs)
wo_tab_file.close()

os.system('cut -f1,3,4,5,6,7,8,9,12,13,14,15,20,21 {} > {}'.format(args.out_readable, args.out_important))
