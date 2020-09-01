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
        '-mr', '--move_range', type=int,
        help = 'search range to find the sequence with lowest GC content, default is 1/4 of the interval'
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


def choose_low_gc_seq(SEQ, probe_length, range_st, range_ed):
    min_GC = 100
    min_idx = None
    min_SEQ = None
    for i in range(range_st, range_ed+1):
        candidate_SEQ = SEQ[i:i+probe_length]
        candidate_GC = gc_percentage(candidate_SEQ)
        if candidate_GC < min_GC:
            min_GC = candidate_GC
            min_idx = i
            min_SEQ = candidate_SEQ
    return (min_idx, min_SEQ)


def design_probe(SEQ, probe_length, probe_number, move_range):
    seq_length = len(SEQ)
    list_probe = []
    if probe_number == 1:
        start_position = round(seq_length - probe_length)/2
        range_st = max(start_position - move_range, 0)
        range_ed = min(start_position + move_range, max(seq_length-probe_length, 0))
        (probe_position, probe_SEQ) = choose_low_gc_seq(SEQ, probe_length, range_st, range_ed)
        list_probe.append((probe_position, probe_SEQ))
    else:
        interval = round((seq_length - probe_length)/(probe_number - 1))
        if move_range == None:
            move_range = max(0, round((interval-probe_length)/4))
        start_position = 0
        for idx in range(probe_number-1):
            range_st = max(0, start_position - move_range) 
            range_ed = start_position + move_range
            probe_position, probe_SEQ = choose_low_gc_seq(SEQ, probe_length, range_st, range_ed) # SEQ[start_position:start_position+probe_length]
            list_probe.append((probe_position, probe_SEQ))
            start_position += interval

        start_position = seq_length - probe_length
        range_st = start_position - move_range
        range_ed = start_position
        (probe_position, probe_SEQ) = choose_low_gc_seq(SEQ, probe_length, range_st, range_ed)
        list_probe.append((probe_position, probe_SEQ))

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
    move_range   = args.move_range
    gc_percentage_threshold = args.gc_percentage_threshold
    fo_probe_csv = args.fo_probe_csv
    
    dict_target = parse_fasta(fn_fasta)
    f_o = open(fo_probe_csv, 'w')
    max_GC = 0
    min_GC = 100
    total_GC = 0
    seq_num  = 0
    gct_num  = 0
    for (name, SEQ) in sorted(dict_target.items()):
        list_probe = design_probe(SEQ, probe_length, probe_number, move_range)
        #print(name)
        #print("seq_length =", len(SEQ))
        try:
            parsed_name = name.split('|')[1]
        except:
            parsed_name = name.split()[0]
        for idx, probe_info in enumerate(list_probe):
            probe_pos = probe_info[0]
            probe_SEQ = probe_info[1]
            GC = gc_percentage(probe_SEQ)
            if GC > gc_percentage_threshold:
                print("WARNING:", parsed_name + '_' + str(idx), "got", round(GC,2), "% GC content.")
                gct_num += 1
            # === output part === #
            f_o.write(parsed_name + '_' + str(idx) + '_pos:' + str(probe_pos) + ',')
            f_o.write(probe_SEQ + ',')
            f_o.write(str(round(GC,2)) + '\n')
            # === output part === #
            if GC > max_GC: max_GC = GC
            if GC < min_GC: min_GC = GC
            total_GC += GC
            seq_num += 1
    f_o.close()
    # report
    print()
    print(seq_num, "probes generated...")
    print("Max GC in probes is", round(max_GC,2), "%. Min GC in probes is", round(min_GC,2), "%. Average GC in probes is", round(total_GC/seq_num,2), "%")
    print("There are", gct_num, "probes exceed threshold", gc_percentage_threshold, "% GC content.")


