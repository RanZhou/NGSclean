#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This code combine counts from feature counts
"""

from sys import argv
import glob
import os

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
# auto

def main():
    script, working_dir, data_out = argv
    out_stream = open(data_out, 'w')
    list_f = glob.glob(os.path.join(working_dir, '*_counts.txt'))
    file_out_list = []
    data_dict = AutoVivification()
    for file_in in sorted(list_f):
        dir2, file_short = os.path.split(file_in)
        if file_short[0].isdigit():
            file_short = 'S'+file_short
        file_short = file_short[0:-4] # remove the suffix
        file_out_list.append(file_short)
        line_count = 0
        with open(file_in, 'r') as in_stream:
            for line in in_stream:
                line_count += 1
                if line_count < 3:
                    continue
                line_record = line.rstrip().split("\t")
                gene_id = line_record[0]
                count = int(round(float(line_record[6])))
                data_dict[gene_id][file_short] = str(count)
    # get the output
    head_line = "gene"+"\t"+"\t".join(file_out_list)
    head_line = head_line.replace('-','_')
    out_stream.write(head_line+"\n")
    for gene_id in sorted(data_dict.keys()):
        out_line = []
        for file_id in file_out_list :
            out_c = '0'
            if file_id in data_dict[gene_id].keys():
                out_c = data_dict[gene_id][file_id]
            out_line.append(out_c)
        out_stream.write(gene_id + "\t" + "\t".join(out_line) + "\n")
    out_stream.close()

## Author : lxue@uga.edu

if __name__ == '__main__':
    main()
