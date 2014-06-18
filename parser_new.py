from Levenshtein import *

from warnings import warn

import swalign

import pickle

import numpy as np

import sklearn

import pdb

# ##################     READS FROM FILE       ########################


def get_reads():
    f = file("/home/alex/another/diploma/parser/ag_1_GGCTAC_filtered.fastq", "r")

    reads_n = []  # read names
    reads_seq = []  # read sequences
    i = 0

    for line in f:
        li = line.split(' ')
        if (i + 1) % 4 == 1:
            reads_n.append(li[0][1:])
        if (i + 1) % 4 == 2:
            reads_seq.append(li[0])
        i += 1

    f.close()
    return reads_n, reads_seq

    # #####################      ALIGNMENTS    #########################


def align():
    o_res = []
    f = open("/home/alex/another/diploma/parser/results_Oases.txt", "r")
    o_res = pickle.load(f)
    f.close()

    t_res = []
    f = open("/home/alex/another/diploma/parser/results_Trinity.txt", "r")
    t_res = pickle.load(f)
    f.close()
    return o_res, t_res

    # ################    READS FOR OASES     #####################


def reads_for_oases():
    f = file("/home/alex/another/diploma/parser/results_oases.sam", "r")

    o_read_name = []  # read name
    o_tr_name = []  # transcript for this read name (if none, transcript is *)
    # o_read_coord = []  # read coordinates in this transcript
    o_read_seq = []

    for line in f:
        columns = line.split('\t')
        o_read_name.append(columns[0])
        o_tr_name.append(columns[2])
        # o_read_coord.append(columns[3])
        o_read_seq.append(columns[9])

    f.close()
    return o_read_name, o_tr_name, o_read_seq

    # ####################       READS FOR TRINITY      #######################


def reads_for_trinity():
    f = file("/home/alex/another/diploma/parser/results_trinity.sam", "r")

    t_read_name = []  # read name
    t_tr_name = []  # transcript for this read name (if none, transcript is *)
    # t_read_coord = []  # read coordinates in this transcript
    t_read_seq = []

    for line in f:
        columns = line.split('\t')
        t_read_name.append(columns[0])
        t_tr_name.append(columns[2])
        # t_read_coord.append(columns[3])
        t_read_seq.append(columns[9])

    f.close()

    return t_read_name, t_tr_name, t_read_seq

    # ########################   ASSEMBLIES READING   ###########################


def assemblies_reading():
    f = open("/home/alex/another/diploma/parser/ref_for_reads.fasta", "r")
    ref = []

    tmp = ""
    while 1:
        l = f.readline()
        if not l:
            ref.append(tmp)
            tmp = ""
            break
        if l[0] != '>':
            tmp += l[:-1]

        else:
            ref.append(tmp)
            tmp = ""

    f.close()

    f = open("/home/alex/another/diploma/parser/Oases.fasta", "r")
    oases = []
    o_name_index = []

    tmp = ""
    while 1:
        l = f.readline()
        if not l:
            oases.append(tmp)
            tmp = ""
            break
        if l[0] != '>':
            tmp += l[:-1]
        else:
            oases.append(tmp)
            tmp = ""
            massive = l.split(' ')
            o_name_index.append(massive[0][1:-1])

    f.close()

    f = open("/home/alex/another/diploma/parser/Trinity.fasta", "r")
    trinity = []
    t_name_index = []

    tmp = ""
    while 1:
        l = f.readline()
        if not l:
            trinity.append(tmp)
            tmp = ""
            break
        if l[0] != '>':
            tmp += l[:-1]
        else:
            trinity.append(tmp)
            tmp = ""
            massive = l.split(' ')
            t_name_index.append(massive[0][1:])

    f.close()
    return ref, oases, o_name_index, trinity, t_name_index


# #####################      DISTANCES     #######################


def dist_read():
    o_dist = np.zeros([len(ref), len(oases)])
    f = open("/home/alex/another/diploma/parser/Levenshtein_for_Oases.txt", "r")
    o_dist = pickle.load(f)
    f.close()

    t_dist = np.zeros([len(ref), len(trinity)])
    f = open("/home/alex/another/diploma/parser/Levenshtein_for_Trinity.txt", "r")
    t_dist = pickle.load(f)
    f.close()

    o_pairs = np.array(range(len(ref)))
    f = open("/home/alex/another/diploma/parser/Similar_transkripts_Oases.txt", "r")
    o_pairs = pickle.load(f)
    f.close()

    t_pairs = np.array(range(len(ref)))
    f = open("/home/alex/another/diploma/parser/Similar_transkripts_Trinity.txt", "r")
    t_pairs = pickle.load(f)
    f.close()

    return o_dist, t_dist, o_pairs, t_pairs


def make_dict_r_t(o_read_name, o_tr_name):
    o_dict_r_t = dict()

    j = 0
    for i in o_read_name:
        o_dict_r_t[i] = o_tr_name[j]
        j += 1

    return o_dict_r_t


def make_dict_t_r(o_name_index, o_tr_name, o_read_name):
    o_dict_t_r = dict()

    for i in o_name_index:
        o_dict_t_r[i] = []

    j = 0
    t = 0

    for i in o_tr_name:
        if o_dict_t_r.has_key(i):
            o_dict_t_r[i].append(o_read_name[j])
        j += 1

    return o_dict_t_r


def make_dict(t_name_index):
    t_dict = dict()

    j = 0
    for i in t_name_index:
        t_dict[i] = j
        j += 1

    return t_dict


def reads_for_class(o_res, o_name_index, o_dict_t_r, t_dict_r_t, t_dict, t_pairs, t_res, reads_dict):
    class_trinity_t = []
    i = 0

    top_bound = 0.9
    bottom_bound = 0.5

    for aln in o_res:
        if i > len(o_name_index):
            break

        print "stage"

        if aln.identity > top_bound:

            transcript_name = o_name_index[i]
            # transcript_seq = oases[i + 1]

            reads_for_tr = []
            if o_dict_t_r.has_key(transcript_name):
                reads_for_tr = o_dict_t_r[transcript_name]
            else:
                pdb.set_trace()
            for read in reads_for_tr:
                if t_dict_r_t[read] == '*':
                    continue
                tr_for_read = t_dict_r_t[read]
                #                pdb.set_trace()
                l = t_dict[tr_for_read]

                m = 0
                u = 0
                for pair in t_pairs:
                    if pair == l:
                        u = 1
                        break
                    m = m + 1

                if u == 1 and t_res[m].identity > bottom_bound:
                    break
                else:
                    class_trinity_t.append(reads_dict[read])
                    # print "succsess!"
                    break

        i += 1
    return class_trinity_t


'''
def making_dict (o_read_name, ):
    o_dict_r_t = dict()

    j = 0
    for i in o_read_name:
        o_dict_r_t[i] = o_tr_name[j]
        j += 1

    t_dict_r_t = dict()

    j = 0
    for i in t_read_name:
        t_dict_r_t[i] = t_tr_name[j]
        j += 1

    o_dict_t_r = dict()

    for i in o_name_index:
        o_dict_t_r[i] = []

    j = 0
    t = 0

    for i in o_tr_name:
        if o_dict_t_r.has_key(i):
            o_dict_t_r[i].append(o_read_name[j])
        j += 1

    t_dict_t_r = dict()

    for i in trinity:
        t_dict_t_r[i] = []

    j = 0

    for i in t_tr_name:
        if t_dict_t_r.has_key(i):
            t_dict_t_r[i].append(t_read_name[j])
        j += 1

    t_dict = dict()

    j = 0
    for i in t_name_index:
        t_dict[i] = j
        j += 1

    return o_dict_r_t, o_dict_t_r, t_dict_r_t, t_dict_t_r, t_dict
'''
if __name__ == '__main__':
    (reads_n, reads_seq) = get_reads()
    if not reads_n or not reads_seq:
        print ("reads have been read unsuccessfully")

    (o_res, t_res) = align()
    if not o_res or not t_res:
        print ("align have been read unsuccessfully")

    (o_read_name, o_tr_name, o_read_seq) = reads_for_oases()
    if not o_read_name or not o_tr_name or not o_read_seq:
        print ("reads for oases have been read unsuccessfully")

    (t_read_name, t_tr_name, t_read_seq) = reads_for_trinity()
    if not t_read_name or not t_tr_name or not t_read_seq:
        print ("reads for trinity have been read unsuccessfully")

    (ref, oases, o_name_index, trinity, t_name_index) = assemblies_reading()
    if not ref or not oases or not trinity or not o_name_index or not t_name_index:
        print ("assemblies have been read unsuccessfully")

    (o_dist, t_dist, o_pairs, t_pairs) = dist_read()
    if t_dist is None or o_dist is None or o_pairs is None or t_pairs is None:
        print ("unsuccessful distance reading")

    o_dict_r_t = make_dict_r_t(o_read_name, o_tr_name)
    t_dict_r_t = make_dict_r_t(t_read_name, t_tr_name)
    o_dict_t_r = make_dict_t_r(o_name_index, o_tr_name, o_read_name)
    t_dict_t_r = make_dict_t_r(t_name_index, t_tr_name, t_read_name)
    o_dict = make_dict(o_name_index)
    t_dict = make_dict(t_name_index)
    reads_dict = make_dict_r_t(reads_n, reads_seq)

    # #########################  READS FOR CLASS TRINITY   #############################

    '''
    class_trinity = []

    i = 0

    top_bound = 0.9
    bottom_bound = 0.5

    for aln in o_res:
        print "stage"

        if aln.identity > top_bound:

            transcript_name = o_name_index[i]
            transcript_seq = oases[i + 1]

            reads_for_tr = []
            if o_dict_t_r.has_key(transcript_name):
                reads_for_tr = o_dict_t_r[transcript_name]
            else:
                pdb.set_trace()
            for read in reads_for_tr:
                if t_dict_r_t[read] == '*':
                    continue
                tr_for_read = t_dict_r_t[read]
#                pdb.set_trace()
                l = t_dict[tr_for_read]

                m = 0
                u = 0
                for pair in t_pairs:
                    if pair == l:
                        u = 1
                        break
                    m = m + 1

                if u == 1 and t_res[m].identity > bottom_bound:
                    break
                else:
                    class_trinity.append(read)
                    # print "succsess!"
                    break

        i += 1
    '''
    '''

    print "o_res[0]"
    print o_res[0]
    print "####################################################"
    print "o_name_index[1]"
    print o_name_index[1]
    print "####################################################"
    print "oases[2]"
    print oases[2]
    print "####################################################"
    print "o_tr_name[0]"
    print o_tr_name[0]
    print "####################################################"
    print "o_read_name[1]"
    print o_read_name[1]
    print "####################################################"
    print "o_read_seq[1]"
    print o_read_seq[1]
    print "####################################################"
    print "t_read_name[1]"
    print t_read_name[1]
    print "####################################################"
    print "t_tr_name[1]"
    print t_tr_name[1]
    print "####################################################"
    print "t_name_index[1]"
    print t_name_index[1]
    print "####################################################"
    print "t_pair[1]"
    print t_pair[1]
    print "####################################################"
    print "t_res[1]"
    print t_res[1]

    '''
    class_trinity_t = reads_for_class(o_res, o_name_index, o_dict_t_r, t_dict_r_t, t_dict, t_pairs, t_res, reads_dict)
    f = file("/home/alex/another/diploma/parser/Class_Trinity.txt", "w")
    pickle.dump(class_trinity_t, f)
    f.close()

    class_trinity_o = reads_for_class(t_res, t_name_index, t_dict_t_r, o_dict_r_t, o_dict, o_pairs, o_res, reads_dict)
    f = file("/home/alex/another/diploma/parser/Class_Oases.txt", "w")
    pickle.dump(class_trinity_o, f)
    f.close()

    print reads_seq[1]
# import cProfile

#cProfile.run('main()')