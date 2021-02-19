import re

fh = open(snakemake.input[0], 'r')
fh_out = open(snakemake.output[0], 'w')

for line in fh:
    if line.startswith('#'):
        continue
    else:
        wo_gaps = re.sub(' +', ' ', line)
        w_tabs = wo_gaps.replace(' ', '\t')
        w_tabs = re.split(r'\t+', w_tabs)
        important = w_tabs[3:9] + w_tabs[12:15] + w_tabs[20:21]
        important = '\t'.join(important)
        fh_out.write('{}\n'.format(important))
fh_out.close()
