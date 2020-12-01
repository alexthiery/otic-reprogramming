#!/opt/conda/bin/python

import re

# create output file to write to
outfile = open("/Users/alex/dev/repos/otic-reprogramming/GalGal6_protein_coding_modified.gtf", 'a')

with open("/Users/alex/dev/repos/otic-reprogramming/output/NF-downstream_analysis/enhancer_analysis/gtf_filter/GalGal6_protein_coding.gtf", 'rt') as gtf:
    for line in gtf:
        # only search lines with gene_id in order to skip header lines
        if 'gene_id' in line:
            if 'gene_name' in line:
                if line[0].isdigit():
                    gene_name = re.sub('.*gene_name "', '', line)
                    gene_name = re.sub('".*', '', gene_name).rstrip()
                    new_line = re.sub('gene_id.*?gene_version', 'gene_id "'+gene_name+'"; gene_version', line)
                    outfile.write('chr' + new_line)
            else:
                if line[0].isdigit():
                    outfile.write('chr' + line)
