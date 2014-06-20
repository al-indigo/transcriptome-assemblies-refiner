from Levenshtein import *

from warnings import warn

import swalign

import pickle

import numpy as np


###################     READS FROM FILE       ########################

f = file("/home/alex/another/diploma/parser/ag_1_GGCTAC_filtered.fastq", "r")

reads_n = []            #read names
reads_seq = []          #read sequences
i = 0

for line in f:
    if (i+1)%4 == 1:
        reads_n.append(line)
    if (i+1)%4 == 3:
        reads_seq.append(line)

f.close()

######################      ALIGNMENTS    #########################

o_res = []
f = open("/home/alex/another/diploma/parser/results_Oases.txt", "r")
o_res = pickle.load(f)
f.close()

t_res = []
f = open("/home/alex/another/diploma/parser/results_Trinity.txt", "r")
t_res = pickle.load(f)
f.close()

#################    READS FOR OASES     #####################

f = file("/home/alex/another/diploma/parser/results_oases.sam", "r")

o_read_name = []          #read name
o_tr_name = []            #transkript for this read name (if none, transkript is *)
o_read_coord = []         #read coordinates in this transcript
o_read_seq = []


for line in f:
    # split the line into a list of column values
    columns = line.split('\t')
    o_read_name.append(columns[0])
    o_tr_name.append(columns[2])
#    o_read_coord.append(columns[3])
    o_read_seq.append(columns[9])

f.close()

#####################       READS FOR TRINITY      #######################
f = file("/home/alex/another/diploma/parser/results_trinity.sam", "r")

t_read_name = []          #read name
t_tr_name = []            #transkript for this read name (if none, transkript is *)
t_read_coord = []         #read coordinates in this transcript
t_read_seq = []

for line in f:
    # split the line into a list of column values
    columns = line.split('\t')
    t_read_name.append(columns[0])
    t_tr_name.append(columns[2])
#    t_read_coord.append(columns[3])
    t_read_seq.append(columns[9])

f.close()
#########################   ASSEMBLIES READING   ###########################

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

######################      DISTANCES     #######################
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



'''
########################    COVERAGE    ###########################


i = 0
j = 0
k = 0

our_transcript = "Locus_1_Transcript_1/1_Confidence_1.000_Length_1805"

for line in o_name_index:
    if line == our_transcript:
        break
    else:
        j = j + 1

print oases[j+1]

score = np.array(range(0, len(oases[j+1])-1))

for elem in o_tr_name:
    if elem == our_transcript:
#        print o_read_coord[k]
#        print int(o_read_coord[k])

        if int(o_read_coord[k])-1>1800:
            print o_read_coord[k]
        for i in range (int(o_read_coord[k])-1, min(int(o_read_coord[k])+len(o_read_seq[i])-1, len(oases[j+1])-1)) :
            score[i] = score[i] + 1
    k = k + 1

for i in range (0, len(oases[j+1])-1):
    print score[i]


#print "len", len(oases[j+1])
    
'''
##########################  READS FOR CLASS TRINITY   #############################

class_trinity = []

i = 0

top_bound = 0.9
bottom_bound = 0.5

for aln in o_res:
    print "stage"

    if aln.identity>top_bound:

        transcript_name = o_name_index[i]
        transcript_seq = oases[i+1]

        j = 1
        for trans in o_tr_name:

            if trans == transcript_name:
                read_name = o_read_name[j]
#                print j
                read_seq = o_read_seq[j]

                k = 0
                for rd in t_read_name:
                    if rd == read_name:
                        tr_an = t_tr_name[k]

                        l = 0
                        for tr_in in t_name_index:
                            if tr_in == tr_an:
                                break
                            else:
                                l = l + 1
#                        print t_pairs.index(l)

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
                            class_trinity.append(read_seq)
#                            print "succsess!"
                            break

                    k = k + 1

            j = j + 1

    i = i + 1
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
f = file("/home/alex/another/diploma/parser/Class_Trinity_2.txt", "w")
pickle.dump (class_trinity, f)
f.close

