#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script trim the adapters of DNA reads, it provide three options to handle the singleton reads.
Author: Liangjiao Xue liangjiao.xue@gmail.com
Copy Right : Tsai Lab, the University of Georgia, cjtsai@uga.edu
             Liangjiao Xue, liangjiao.xue@gmail.com
"""


from optparse import OptionParser
import glob
import os
import gzip
import sys


class ReadRecord:
    def __init__(self,thread):
        self.paired = False
        self.thread = thread
        self.samples = ()
        self.sample2len = dict()
        self.sample2phred = dict()
        self.sample2end_tag = dict()  # 'end', "mid', 'no'


    def read_design(self, design_file):
        with open(design_file, 'r') as design_f:
            lines = [line.rstrip().split("\t") for line in design_f]
            if len(lines[0]) == 3:
                self.paired = True
            elif len(lines[0]) == 2:
                self.paired = False
            else:
                print("Incorrect design format, please double check the file")
            if self.paired :
                if lines[0][2] != 'right':
                    sys.exit("Please check the head of design files")
            self.samples = lines[1:]

    def check_fastq(self):
        for sample_record in self.samples:
            first_fqgz = sample_record[1]
            sample_name = sample_record[0]
            sample_name.replace(' ', '_')
            # get the phred format
            len_max = 0
            qual_min = 1000
            qual_max = 0
            # changes by Jake to have try/finally block instead of for loop. Fixes GzipFile error.
            fqgzIN = gzip.open(first_fqgz, "rt")
            try:
                count = 0
                limit = 200
                qual_str_all = ''
                for line in fqgzIN:
                    count += 1
                    if count == 1:
                        name_line = line.rstrip()
                        tag_name = 'no'
                        if name_line.endswith('/1'):
                            tag_name = 'end'
                        elif ' 1:' in name_line:
                            tag_name = 'mid'
                        self.sample2end_tag[sample_name] = tag_name
                    if count > limit:
                        break
                    if count % 4 == 0:
                        qual_str = line.rstrip()
                        qual_str_all += qual_str
                        if len_max < len(qual_str):
                            len_max = len(qual_str)
                for qual_char in qual_str_all:
                    qual = ord(qual_char)
                    if qual_min > qual:
                        qual_min = qual
                    if qual_max < qual:
                        qual_max = qual
            finally:
                fqgzIN.close()
            if qual_min < 33 or qual_max > 105:
                sys.exit("Quality values corrupt. found [$min; $max] where [33; 104] was expected \n")
            if qual_min >= 59 and qual_max <= 110:
                qual_out = '64'
            elif qual_min >= 33 and qual_max <= 74:
                qual_out = '33'
            else:
                sys.exit("May be new fastq format \n")
            # output to class
            self.sample2len[sample_name] = len_max
            self.sample2phred[sample_name] = qual_out



def main():
    ########################################################################
    # Options
    parser = OptionParser()
    usage = "Usage: %prog [options] arg1 arg2"
    parser.add_option("-d", "--design", dest="design",
                      help="Input design file")
    parser.add_option("-s", "--singleton", dest="singleton", default="keep",
                      help="How to treat singleton reads: keep(default), merge with paired, or discard")
    parser.add_option("-t", "--threadNum", dest="thread",
                      help="Numbers of thread")
    parser.add_option("-o", "--outputDir", dest="output",
                      help="Change the workingDir for read output")
    parser.add_option("-q", "--queue", dest="queue",
                      help="Queue for batch jobs, inter (run interactively) or queue (default)")
    parser.add_option("-n", "--node", dest="node",
                      help="Node of queue")
    parser.add_option("-z", "--run_trimmomatic", dest="trimmomatic",
                      help="path to trimmomatic.jar")
    parser.add_option("-x", "--load_trimmo_module", dest="loadtrimmo",
                      help="module for trimmomatic")
    parser.add_option("-a", "--adaptor", dest="adaptor",
                      help="path to trimmomatic adaptor")


    (options, args) = parser.parse_args()
    # parser options
    option_dict = vars(options)
    design = option_dict['design']
    singleton  = option_dict['singleton']
    thread = option_dict['thread']
    output = option_dict['output']
    queue  = option_dict['queue']
    node   = option_dict['node']
    adaptor= option_dict['adaptor']
    trimmomatic = option_dict['trimmomatic']
    loadtrimmo = option_dict['loadtrimmo']

    #
    if thread is None:
        thread = 1
    if queue is None:
        queue = 'batch'
    if node is None:
        node = 'Intel'
    if trimmomatic is None :
        sys.exit("Please provide the path to trimmomatic.jar")
    if adaptor is None:
        sys.exit("Please provide the path to trimmomatic adaptor")
    if output is None:
        output = os.getcwd()
    python_dir_path = os.path.dirname(os.path.realpath(__file__))
    #print(python_dir_path +"path to python code")
    ###########################################################################
    # read design file
    readrecord = ReadRecord(thread)
    readrecord.read_design(design)
    readrecord.check_fastq()

    # prepare shell script
    dir, file_design_short = os.path.split(design)
    master_file = 'Run_'+file_design_short[0:-4]+'.sh'
    working_dir = output

    master_out = open(master_file, 'w')
    master_out.write("#!/bin/sh\n")
    master_out.write("cd " + working_dir + "\n")

    #############################################################################
    ## Each sample
    count = 0
    for sample_record in readrecord.samples:
        ## get children shell name
        sample_name = sample_record[0]
        sample_name.replace(' ', '_')
        first_fastq = sample_record[1]
        first_fastq.replace(' ', '')
        dir, first_fastq_short = os.path.split(first_fastq)
        prefix = sample_name

        count += 1
        file_shell_str = 'r' + str(count) + '_' + prefix + '.sh'
        if queue == 'inter':
            master_out.write("chmod 750 " + file_shell_str + "\n")
            master_out.write("./" + file_shell_str + "\n")
        else:
            master_out.write("qsub " + file_shell_str + "\n")

        #############################
        ## get log file
        final_log = prefix + "_Final.log.txt"
        log_out_stream = open(final_log,"w")
        log_out_stream.write("Read length:" + str(readrecord.sample2len[sample_name])+"\n")
        log_out_stream.write("Read Qual foramt:" + str(readrecord.sample2phred[sample_name]) + "\n")
        log_out_stream.close()

        ##############################
        ## write childredn shell
        if queue == 'inter':
            file_shell = open(file_shell_str, 'w')
            file_shell.write("cd " + working_dir + "\n\n")
            file_shell.close()
        else:
            write_shell_head(file_shell_str, queue, thread, node, working_dir)

        ##############################
        ## write trmmomatic shell
        all_fastq_file = " ".join(sample_record[1:])
        write_shell_trimmo(file_shell_str, loadtrimmo, trimmomatic, prefix, all_fastq_file,readrecord,sample_name,adaptor)

        ### put the log here

        ###################################################
        ### Summary data
        final_log = prefix + "_Final.log.txt"
        trim_log = prefix + '_trim.log'
        # out read summary first
        log_out_stream = open(final_log, "a")
        write_shell_summary(file_shell_str,final_log,trim_log)

        ###################################################
        ###  prepare fastq output
        paired_boolean = readrecord.paired
        end_tag = readrecord.sample2end_tag[sample_name]
        write_shell_fastq(file_shell_str, paired_boolean, prefix, singleton, end_tag)

        ##################################################
        ### zip and clean
        write_shell_clean(file_shell_str,prefix, paired_boolean,singleton)
# End of main function
################################################################################################################

def write_shell_summary(file_shell_str, final_log, trim_log):
    file_shell = open(file_shell_str, 'a')
    file_shell.write('printf "Trimming : " >>' + final_log + "\n")
    file_shell.write('grep "^Input Read"  ' + trim_log + '>>' + final_log + "\n")
    file_shell.close()


def write_shell_clean(file_shell_str, prefix, paired_boolean, singleton):
    file_shell = open(file_shell_str, 'a')
    remove_list = [prefix + '_trim.log',]
    if paired_boolean:
        # paired
        if singleton == 'discard':
           fastq_temp_list = ( "_trimS_1.fq.gz", "_trimS_1.fq.gz")
           for suffix in fastq_temp_list:
                remove_list.append(prefix + suffix)
        elif singleton == 'merge':
            fastq_temp_list = ("_trimP_1.fq.gz","_trimP_2.fq.gz","_trimS_1.fq.gz", "_trimS_2.fq.gz",
                               "_clean_F1.fq.gz","_clean_F2.fq.gz")
            for suffix in fastq_temp_list:
                remove_list.append(prefix + suffix)
    ## Write shell
    for remove_f in remove_list:
        file_shell.write("rm " + remove_f + "\n")
    file_shell.close()


def write_shell_head(file_shell_str, queue,thread,  node, working_dir):
    # children shell  head
    file_shell = open(file_shell_str, 'w')
    file_shell.write("#!/bin/sh\n" +
                "#PBS -N j-clean" + file_shell_str[0:-3] + "\n" +
                "#PBS -q " + queue + "\n" +
                "#PBS -l nodes=1:ppn=" + str(thread) + ':' + node + "\n" +
                "#PBS -l walltime=48:00:00\n" +
                "#PBS -l mem=20gb\n" +
                "cd " + working_dir + "\n\n")
    file_shell.close()

def write_shell_trimmo(file_shell_str, loadtrimmo, trimmomatic, prefix, all_fastq_file,readrecord,sample_name,adaptor):
    file_shell = open(file_shell_str, 'a')
    paired_boolean = readrecord.paired
    phred_tag = ''
    if readrecord.sample2phred[sample_name] == '64':
        phred_tag = '-phred64'
    else:
        phred_tag = '-phred33'
    paired_boolean = readrecord.paired
    thread = readrecord.thread
    ### write
    if loadtrimmo is None:
        print("trimmomatic module is not loaded")
    else:
        file_shell.write("module load " + loadtrimmo + " \n")
    pe_tag = 'SE'
    trim_out = ''
    trim_log = prefix + '_trim.log'
    if paired_boolean:
        # Paired end
        pe_tag = 'PE'
        r1_pair = prefix + "_trimP_1.fq.gz"
        r1_unpair = prefix + "_trimS_1.fq.gz"
        r2_pair = prefix + "_trimP_2.fq.gz"
        r2_unpair = prefix + "_trimS_2.fq.gz"
        trim_out_list = (r1_pair, r1_unpair, r2_pair, r2_unpair)
        trim_out = " ".join(trim_out_list)
    else:        # single end
        pe_tag = 'SE'
        trim_out = prefix + '_trim_1.fq.gz'
    file_shell.write("time java -jar " + trimmomatic + ' ' + pe_tag + " -threads " + str(thread) + " " +
                     phred_tag + " \\\n" +
                     " " + all_fastq_file + " \\\n" +
                     " " + trim_out + " \\\n" +
                     " ILLUMINACLIP:" + adaptor + ':2:30:10' +
                     ' LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 &>' + trim_log + "\n\n")
    file_shell.close()



def write_shell_fastq(file_shell_str, paired_boolean, prefix, singleton, end_tag):
    file_shell = open(file_shell_str, 'a')
    if not paired_boolean:
        # single send
        trim_out = prefix + '_trim_1.fq.gz'
        final_out1 = prefix + "_clean_1.fq.gz"
        file_shell.write("mv " + trim_out +" "+ final_out1 + "\n")
    else:
        # Paired end
        final_out1 = prefix + "_clean_1.fq.gz"
        final_out2 = prefix + "_clean_2.fq.gz"
        r1_pair = prefix + "_trimP_1.fq.gz"
        r2_pair = prefix + "_trimP_2.fq.gz"
        if singleton == 'merge':
            # regenerate records
            r1_unpair = prefix + "_trimS_1.fq.gz"
            r2_unpair = prefix + "_trimS_2.fq.gz"
            fake_fq2_from_fq1 = prefix + "_clean_F2.fq"
            if end_tag == 'no' :
                file_shell.write(''' zcat ''' + r1_unpair + ' | ' + " \\\n" +
                        ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                             fake_fq2_from_fq1 + "\n")
            elif end_tag == 'end' :
                file_shell.write('zcat '+r1_unpair + '|' + ''' sed '1~4 s|/1$|/2|g'  ''' + ' | ' + " \\\n" +
                         ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                             fake_fq2_from_fq1 + "\n")
            elif end_tag == 'mid' :
                file_shell.write('zcat '+r1_unpair + '|' + ''' sed '1~4 s|\s1\:| 2\:|g' ''' + ' | ' + " \\\n" +
                        ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                            fake_fq2_from_fq1 + "\n")
            else :
                print("New fastq format\n")
            # fake 2
            fake_fq1_from_fq2 = prefix + "_clean_F1.fq"
            if end_tag == 'no':
                file_shell.write(''' zcat ''' + r2_unpair + ' | ' + " \\\n" +
                                 ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                                 fake_fq1_from_fq2 + "\n")
            elif end_tag == 'end':
                file_shell.write('zcat '+r2_unpair+ '|'+''' sed '1~4 s|/2$|/1|g'  '''   + ' | ' + " \\\n" +
                                 ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                                 fake_fq1_from_fq2 + "\n")
            elif end_tag == 'mid':
                file_shell.write('zcat '+r2_unpair+ '|'+''' sed '1~4 s|\s2\:| 1\:|g'  '''  + ' | ' + " \\\n" +
                                 ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                                 fake_fq1_from_fq2 + "\n")
            else:
                print("New fastq format\n")
            # merge the files
            file_shell.write('gzip ' + fake_fq2_from_fq1 + "\n")
            file_shell.write('gzip ' + fake_fq1_from_fq2 + "\n")
            file_shell.write('cat ' + r1_pair + " " + r1_unpair + " " + fake_fq1_from_fq2 + '.gz' + " >" + final_out1 + "\n")
            file_shell.write('cat ' + r2_pair + " " + fake_fq2_from_fq1 + '.gz' + " " + r2_unpair + " >" + final_out2 + "\n")
        else :
            file_shell.write('mv  ' + r1_pair + " " + final_out1 + "\n")
            file_shell.write('mv  ' + r2_pair + " " + final_out2 + "\n")
    file_shell.close()

#######################################################
#  Main function
#######################################################

if __name__ == '__main__':
    main()
