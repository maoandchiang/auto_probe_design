import argparse
import pickle
import os
import numpy as np
#from utils import get_reverse_complement
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-fa', '--fn_fasta',
        help = 'input target sequence fasta file'
    )
    parser.add_argument(
        '-pl', '--probe_length', type=int, default=60,
        help = 'design probe length'
    )
    parser.add_argument(
        '-pn', '--probe_number', type=int, default=1,
        help = 'design probe number'
    )
    parser.add_argument(
        '-gct', '--gc_percentage_threshold', type=int, default=60,
        help = 'the gc content percentage threshold'
    )
    
    parser.add_argument(
        '-fo', '--fo_probe_csv',
        help = 'output designed probe csv file'
    )
    args = parser.parse_args()
    return args

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def parse_fasta(fn_fasta):
    '''parse the fasta file into a dictionary'''
    # dict_name_SEQ {}
    #  - keys: seq_name
    #  - values: seq_SEQ
    dict_name_SEQ = {}
    with open(fn_fasta, 'r') as f_f:
        seq_name = ""
        seq_SEQ = ""
        for line in f_f:
            if line[0] == '>':
                if seq_name != "":
                    if dict_name_SEQ.get(seq_name):
                        print("WARNING! Duplicate sequence name:", seq_name)
                        new_name = assign_new_name_basic(seq_name, dict_name_SEQ)
                        dict_name_SEQ[new_name] = seq_SEQ
                    else:
                        dict_name_SEQ[seq_name] = seq_SEQ
                seq_name = line.strip()[1:]
                seq_SEQ = ""
            else:
                seq_SEQ += line.strip()
        if dict_name_SEQ.get(seq_name):
            new_name = assign_new_name_basic(seq_name, dict_name_SEQ)
            dict_name_SEQ[new_name] = seq_SEQ
            print("WARNING! Duplicate sequence name:", seq_name)
        else:
            dict_name_SEQ[seq_name] = seq_SEQ
    return dict_name_SEQ


def design_probe(SEQ, probe_length, probe_number):
    seq_length = len(SEQ)
    list_probe = []
    if probe_number == 1:
        start_position = max( (seq_length - probe_length)/2 , 0)
        probe_SEQ = SEQ[start_position:start_position+min(seq_length,probe_length)]
        list_probe.append((start_position, probe_SEQ))
    else:
        interval = int((seq_length - probe_length)/(probe_number - 1))
        start_position = 0
        for idx in range(probe_number-1):
            probe_SEQ = SEQ[start_position:start_position+probe_length]
            list_probe.append((start_position, probe_SEQ))
            start_position += interval
        probe_SEQ = SEQ[-probe_length:]
        start_position = seq_length - probe_length
        list_probe.append((start_position, probe_SEQ))

    return list_probe


def gc_percentage(SEQ):
    SEQ = SEQ.upper()
    percentage = (SEQ.count('C') + SEQ.count('G'))*100/len(SEQ)
    return percentage



if __name__ == "__main__":
    args = parse_args()
    fn_fasta = args.fn_fasta
    probe_length = args.probe_length
    probe_number = args.probe_number
    gc_percentage_threshold = args.gc_percentage_threshold
    fo_probe_csv = args.fo_probe_csv
    
    dict_target = parse_fasta(fn_fasta)
    f_o = open(fo_probe_csv, 'w')
    for (name, SEQ) in sorted(dict_target.items()):
        list_probe = design_probe(SEQ, probe_length, probe_number)
        #print(name)
        #print("seq_length =", len(SEQ))
        parsed_name = name.split('|')[1]
        for idx, probe_info in enumerate(list_probe):
            probe_pos = probe_info[0]
            probe_SEQ = probe_info[1]
            f_o.write(parsed_name + '_' + str(idx) + '_pos:' + str(probe_pos) + ',')
            f_o.write(probe_SEQ + ',')
            f_o.write(str(gc_percentage(probe_SEQ)) + '\n')
    f_o.close()




