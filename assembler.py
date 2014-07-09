import pickle
import nwalign
from tr_parser import *
import time
import Levenshtein
import pdb
import multiprocessing
import glob
import munkres
import numpy as np
from numpy import median
from sklearn.utils.linear_assignment_ import linear_assignment
from sklearn.utils.validation import check_arrays


def get_cost_matrix(first_name_list, second_name_list, first_to_second_similarity_tuples):

    offset = -min([x[3] for x in first_to_second_similarity_tuples]) + 1

    data_dict = {first_name: dict() for first_name in first_name_list}

    for first_transcript_name, \
        second_transcript_name, \
        alignment, \
        score, \
        score_self_first, \
        score_self_second, \
        distance, \
        len_first, \
        len_second in first_to_second_similarity_tuples:

        data_dict[first_transcript_name][second_transcript_name] = \
            1.0 - (score + offset)/(float(max(score_self_first, score_self_second) + offset))

    cost_matrix = []

    for first_tr_name in first_name_list:
        row_dict = data_dict[first_tr_name]

        cost_matrix.append([
            row_dict[second_tr_name] for second_tr_name in second_name_list])

    return cost_matrix


def get_classification_results(filename):
    results = dict()
    with open(filename, 'r') as f:
        try:
            while f:
                results.update(pickle.load(f))
        except EOFError:
            #print len(results)
            f.close()
    return results


def get_pickled_list(filename):
    results = []
    with open(filename, 'r') as f:
        try:
            while f:
                results.extend(pickle.load(f))
        except EOFError:
            #print ('Read: ', len(results))
            f.close()

    return results


def align_thread(first_name_index, second_name_index, first_transcript_by_name, second_transcript_by_name, fileout_pre):
    print multiprocessing.current_process()._identity[0]
    #scoring = swalign.NucleotideScoringMatrix(2, -1)
    #sw = swalign.LocalAlignment(scoring)
    start = time.time()
    iteration = 0
    for first_transcript_name in first_name_index:
        each_with_each_alignment = []
        for second_transcript_name in second_name_index:
            alignment = nwalign.global_align(first_transcript_by_name[first_transcript_name],
                                             second_transcript_by_name[second_transcript_name],
                                             matrix='BLOSUM50.txt')
            score = nwalign.score_alignment(alignment[0],
                                            alignment[1],
                                            gap_open=-5,
                                            gap_extend=-2,
                                            matrix='BLOSUM50.txt')
            score_self_second = nwalign.score_alignment(alignment[0],
                                                         alignment[0],
                                                         gap_open=-5,
                                                         gap_extend=-2,
                                                         matrix='BLOSUM50.txt')
            score_self_first = nwalign.score_alignment(alignment[1],
                                                       alignment[1],
                                                       gap_open=-5,
                                                       gap_extend=-2,
                                                       matrix='BLOSUM50.txt')
            distance = Levenshtein.distance(first_transcript_by_name[first_transcript_name], second_transcript_by_name[second_transcript_name])
            each_with_each_alignment.append((first_transcript_name,
                                             second_transcript_name,
                                             alignment,
                                             score,
                                             score_self_first,
                                             score_self_second,
                                             distance,
                                             len(first_transcript_by_name[first_transcript_name]),
                                             len(second_transcript_by_name[second_transcript_name])))
        iteration += 1
        print (iteration, len(first_name_index), int(time.time() - start))
        out = open(fileout_pre + str(multiprocessing.current_process()._identity[0]), 'wb')
        pickle.dump(each_with_each_alignment, out)
        out.close()


def get_align():
#    f = open('data/oases_fit_transcripts', 'r')
#    oases_fit_transcripts = pickle.load(f)
#    f.close()
#    f = open('data/oases_not_fit_transcripts', 'r')
#    oases_not_fit_transcripts = pickle.load(f)
#    f.close()
#    f = open('data/trinity_fit_transcripts', 'r')
#    trinity_fit_transcripts = pickle.load(f)
#    f.close()
#    f = open('data/trinity_not_fit_transcripts', 'r')
#    trinity_not_fit_transcripts = pickle.load(f)
#    f.close()

    (ref_transcripts, oases_transcripts, oases_name_index, trinity_transcripts, trinity_name_index) = get_assemblies("data/ref_for_reads.fasta",
                                                                                                                     "data/Oases.fasta",
                                                                                                                     "data/Trinity.fasta")
    # it's a hack
    oases_transcripts.pop(0)
    trinity_transcripts.pop(0)
    ref_transcripts.pop(0)

    oases_tr_by_name = dict(zip(oases_name_index, oases_transcripts))
    trinity_tr_by_name = dict(zip(trinity_name_index, trinity_transcripts))

    ref_name_index = range(0, len(ref_transcripts)-1)
    ref_tr_by_name = dict(zip(ref_name_index, ref_transcripts))

    lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]

    tr_parts = lol(trinity_name_index, 21)

    if False:
        process_list = []
        for tr_part in tr_parts:
            p = multiprocessing.Process(target=align_thread, args=(tr_part,
                                                                   oases_name_index,
                                                                   trinity_tr_by_name,
                                                                   oases_tr_by_name,
                                                                   'results/trinity_to_oases_each_align'))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        ref_parts = lol(ref_name_index, 21)

        process_list = []
        for ref_part in ref_parts:
            p = multiprocessing.Process(target=align_thread, args=(ref_part,
                                                                   oases_name_index,
                                                                   ref_tr_by_name,
                                                                   oases_tr_by_name,
                                                                   'results/ref_to_oases_each_align'))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        process_list = []
        for ref_part in ref_parts:
            p = multiprocessing.Process(target=align_thread, args=(ref_part,
                                                                   trinity_name_index,
                                                                   ref_tr_by_name,
                                                                   trinity_tr_by_name,
                                                                   'results/ref_to_trinity_each_align'))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

    trinity_to_oases = []
    ref_to_oases = []
    ref_to_trinity = []
    for filename in glob.glob('./results/sim/trinity_to_oases_each_align*'):
        trinity_to_oases.extend(get_pickled_list(filename))
    for filename in glob.glob('./results/sim/ref_to_oases_each_align*'):
        ref_to_oases.extend(get_pickled_list(filename))
    for filename in glob.glob('./results/sim/ref_to_trinity_each_align*'):
        ref_to_trinity.extend(get_pickled_list(filename))

    print (len(trinity_to_oases), len(ref_to_oases), len(ref_to_trinity))

    offset_ref_oases = -min([x[3] for x in ref_to_oases]) + 1
    offset_ref_trinity = -min([x[3] for x in ref_to_trinity]) + 1

    print (offset_ref_oases, offset_ref_trinity)

    median_score_ref_oases = median([x[3] + offset_ref_oases for x in ref_to_oases])
    median_score_ref_trinity = median([x[3] + offset_ref_trinity for x in ref_to_trinity])

    print(median_score_ref_oases, median_score_ref_trinity)

    median_best_ref_oases1 = median([x[4] + offset_ref_oases for x in ref_to_oases])
    median_best_ref_oases2 = median([x[5] + offset_ref_oases for x in ref_to_oases])
    median_best_ref_oases3 = median([max(x[4], x[5]) + offset_ref_oases for x in ref_to_oases])

    median_best_ref_trinity1 = median([x[4] + offset_ref_trinity for x in ref_to_trinity])
    median_best_ref_trinity2 = median([x[5] + offset_ref_trinity for x in ref_to_trinity])
    median_best_ref_trinity3 = median([max(x[4], x[5]) + offset_ref_trinity for x in ref_to_trinity])

    print (median_best_ref_oases1, median_best_ref_oases2, median_best_ref_oases3)
    print (median_best_ref_trinity1, median_best_ref_trinity2, median_best_ref_trinity3)

    start = time.time()
    print ('Getting cost matrix', time.time() - start)
    ref_to_oases_cost_matrix = get_cost_matrix(ref_name_index, oases_name_index, ref_to_oases)
    print ('Cost matrix ready', time.time() - start)
    indexes = linear_assignment(ref_to_oases_cost_matrix)
    print ('Linear assignment ready', time.time() - start)

    ref_to_oases_correspondence_list = [(ref_name_index[index[0]],
                                         oases_name_index[index[1]],
                                         1.0 - ref_to_oases_cost_matrix[index[0]][index[1]]) for index in indexes]
    f = open('results/ref_to_oases_correspondence_list', 'w')
    pickle.dump(ref_to_oases_correspondence_list, f)
    f.close()
    print ('Linear assignments written', time.time() - start)
    print ('Ref to Oases, median:',
           median([x[2] for x in ref_to_oases_correspondence_list]),
           'Max:',
           max([x[2] for x in ref_to_oases_correspondence_list]),
           'Min:',
           min([x[2] for x in ref_to_oases_correspondence_list]))

    ref_to_trinity_cost_matrix = get_cost_matrix(ref_name_index, trinity_name_index, ref_to_trinity)
    indexes = linear_assignment(ref_to_trinity_cost_matrix)
    ref_to_trinity_correspondence_list = [(ref_name_index[index[0]],
                                           trinity_name_index[index[1]],
                                           1.0 - ref_to_trinity_cost_matrix[index[0]][index[1]]) for index in indexes]
    f = open('results/ref_to_trinity_correspondence_list', 'w')
    pickle.dump(ref_to_trinity_correspondence_list, f)
    f.close()

    print ('Ref to Trinity, median:',
           median([x[2] for x in ref_to_trinity_correspondence_list]),
           'Max:',
           max([x[2] for x in ref_to_trinity_correspondence_list]),
           'Min:',
           min([x[2] for x in ref_to_trinity_correspondence_list]))


    trinity_to_oases_cost_matrix = get_cost_matrix(trinity_name_index, oases_name_index, trinity_to_oases)
    indexes = linear_assignment(trinity_to_oases_cost_matrix)
    trinity_to_oases_correspondence_list = [(trinity_name_index[index[0]],
                                             oases_name_index[index[1]],
                                             1.0 - trinity_to_oases_cost_matrix[index[0]][index[1]]) for index in indexes]
    f = open('results/trinity_to_oases_correspondence_list', 'w')
    pickle.dump(trinity_to_oases_correspondence_list, f)
    f.close()

    print ('Trinity to Oases, median:',
           median([x[2] for x in trinity_to_oases_correspondence_list]),
           'Max:',
           max([x[2] for x in trinity_to_oases_correspondence_list]),
           'Min:',
           min([x[2] for x in trinity_to_oases_correspondence_list]))

    pdb.set_trace()


get_align()

