# -*- coding: utf-8 -*-
import csv, pdb
import pprint
if True:
    def check_validation_results(assembler):
        missing = 0
        existing = 0
        good = 0
        good_list = []
        kernels = [{'svc': 'rbf'}, {'svc': 'poly'}, {'svc': 'sigmoid'}, {'rfc': 2}, {'rfc': 4}, {'rfc': 10}, {'rfc': 100}, {'rfc': 200}, {'rfc': 1000}]
        for bottom_bound in [0.3, 0.4, 0.5, 0.6]:
            for top_bound in [0.99, 0.95, 0.9]:
                for ngram in [4, 5, 6, 7, 8]:
                    for kernel in kernels:
                        matched_must_fit = False
                        matched_must_not_fit = False
                        for fit in ['must-fit', 'must-not-fit']:
                            filename = 'validation/' + assembler + "_classified-" \
                                       + str(ngram) + "-" + kernel.keys()[0] + '-' + str(kernel[kernel.keys()[0]]) + '--' \
                                       + str(bottom_bound) + '-' + str(top_bound) + '--validating-' + fit + "-SUMMARY"
                            try:
                                f = open(filename, 'r')
                                cr = csv.reader(f, delimiter='\t')
                                for row in cr:
                                    if fit == 'must-fit':
                                        #pdb.set_trace()
                                        if 8*int(row[1]) > 9*int(row[2]):
                                            matched_must_fit = True
                                            class_one = [row[1], row[2]]
                                    if fit == 'must-not-fit':
                                        if 8*int(row[2]) > 9*int(row[1]):
                                            matched_must_not_fit = True
                                            class_two = [row[1], row[2]]
                                f.close()
                                existing += 1
                            except IOError:
                                missing += 1
                        if matched_must_fit and matched_must_not_fit:
                            good += 1
                            good_list.append({'bottom': bottom_bound, 'top': top_bound, 'n': ngram,
                                              'kernel': kernel, 'class_one': class_one, 'class_two': class_two})
        print (existing, missing, good)
        pp = pprint.PrettyPrinter()
        pp.pprint(good_list)
        return good_list


    good_list_oases = check_validation_results('oases')
    print "$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    good_list_trinity = check_validation_results('trinity')
    print "$$$$$$$$$$$$$$$$$$$$$$$$$$$"

    for oases_item in good_list_oases:
        for trinity_item in good_list_trinity:
            if oases_item['bottom'] == trinity_item['bottom'] \
                    and oases_item['top'] == trinity_item['top'] \
                    and oases_item['n'] == trinity_item['n']\
                    and oases_item['kernel'].values()[0] == trinity_item['kernel'].values()[0]:
                print (oases_item, trinity_item)

    #intersect = set(good_list_oases) & set(good_list_trinity)
    exit(1)

from classifier import *
from tr_parser import *
from metrics import *
from ngram import NGram

top_bound = 0.9
bottom_bound = 0.5

(ref, oases_reads, oases_name_index, trinity_reads, trinity_name_index) = get_assemblies("data/ref_for_reads.fasta",
                                                                                         "data/Oases.fasta",
                                                                                         "data/Trinity.fasta")

(reads_names, reads_seq) = get_reads("data/ag_1_GGCTAC_filtered.fastq")
if not reads_names or not reads_seq:
    print ("reads have been read unsuccessfully")

(oases_alignment_data, trinity_alignment_data) = get_alignment_data("data/results_Oases.txt",
                                                                    "data/results_Trinity.txt")
if not oases_alignment_data or not trinity_alignment_data:
    print ("align have been read unsuccessfully")

(oases_reads_names, oases_transcripts_names, oases_reads_seq) = get_reads_for_assembler("data/results_oases.sam")
if not oases_reads_names or not oases_transcripts_names or not oases_reads_seq:
    print ("reads for oases have been read unsuccessfully")

(trinity_reads_names, trinity_transcripts_names, trinity_reads_seq) = get_reads_for_assembler("data/results_trinity.sam")
if not trinity_reads_names or not trinity_transcripts_names or not trinity_reads_seq:
    print ("reads for trinity have been read unsuccessfully")

(ref, oases_reads, oases_name_index, trinity_reads, trinity_name_index) = get_assemblies("data/ref_for_reads.fasta",
                                                                                         "data/Oases.fasta",
                                                                                         "data/Trinity.fasta")
if not ref or not oases_reads or not trinity_reads or not oases_name_index or not trinity_name_index:
    print ("assemblies have been read unsuccessfully")

(oases_distances_pairs, trinity_distances_pairs) = get_distances("data/Similar_transkripts_Oases.txt",
                                                                 "data/Similar_transkripts_Trinity.txt")
if oases_distances_pairs is None or trinity_distances_pairs is None:
    print ("unsuccessful distance reading")

reference_reads_to_transcripts = make_index_reads_to_transcripts(reads_names,
                                                                 reads_seq)
oases_reads_to_transcripts = make_index_reads_to_transcripts(oases_reads_names,
                                                             oases_transcripts_names)
trinity_reads_to_transcripts = make_index_reads_to_transcripts(trinity_reads_names,
                                                               trinity_transcripts_names)
oases_transcripts_to_reads = make_index_transcripts_to_reads(oases_name_index,
                                                             oases_transcripts_names,
                                                             oases_reads_names)
trinity_transcripts_to_reads = make_index_transcripts_to_reads(trinity_name_index,
                                                               trinity_transcripts_names,
                                                               trinity_reads_names)
oases_index_by_name = make_index_by_name(oases_name_index)
trinity_index_by_name = make_index_by_name(trinity_name_index)

#TODO: check if it's really class trinity
class_good_trinity = reads_for_class(oases_alignment_data,
                                     oases_name_index,
                                     oases_transcripts_to_reads,
                                     trinity_reads_to_transcripts,
                                     trinity_index_by_name,
                                     trinity_distances_pairs,
                                     trinity_alignment_data,
                                     reference_reads_to_transcripts,
                                     top_bound,
                                     bottom_bound)
f = file("data/Class_Trinity.txt", "w")
pickle.dump(class_good_trinity, f)
f.close()

class_good_oases = reads_for_class(trinity_alignment_data,
                                   trinity_name_index,
                                   trinity_transcripts_to_reads,
                                   oases_reads_to_transcripts,
                                   oases_index_by_name,
                                   oases_distances_pairs,
                                   oases_alignment_data,
                                   reference_reads_to_transcripts,
                                   top_bound,
                                   bottom_bound)
f = file("data/Class_Oases.txt", "w")
pickle.dump(class_good_oases, f)
f.close()


def generate_dictionary(alphabet, length):
    c = [[]]
    dictionary = []
    for i in range(length):
        c = [[x]+y for x in alphabet for y in c]
    for sym_list in c:
        dictionary.append(''.join(sym_list))
    return dictionary


def init_grams_dict(n, alphabet):
    dictionary = generate_dictionary(alphabet, n)
    full_dict = dict()
    for word in dictionary:
        full_dict[word] = float(0)
    return full_dict


def get_distr(strlist, n_len):
    alphabet = ['A', 'C', 'G', 'T', 'N']
    n = NGram(N=n_len, pad_len=0)
    all_ngrams = 0
    grams = init_grams_dict(n_len, alphabet)
    for item in strlist:
        if item == '':
            continue
        ngram_list = list(n._split(item))
        for ng in ngram_list:
            if ng in grams:
                grams[ng] += float(1)
                all_ngrams += 1
    for item in grams.keys():
        grams[item] /= all_ngrams
    return grams

pp = pprint.PrettyPrinter()
#pp.pprint(grams)
#print len(grams)

#pdb.set_trace()
import plotly.plotly as py
from plotly.graph_objs import *
from scipy.stats import wilcoxon

for n in [2, 3, 4, 5, 6, 7, 8]:
    if False:
        reference_distr = get_distr(ref, n)
        full_reads_distr = get_distr(reads_seq, n)
    oases_distr = sorted(get_distr(class_good_oases, n).items())
    oases_bad_distr = sorted(get_distr(list(set([item for item in reads_seq if item not in class_good_oases])), n).items())
    trinity_distr = sorted(get_distr(class_good_trinity, n).items())
    trinity_bad_distr = sorted(get_distr(list(set([item for item in reads_seq if item not in class_good_trinity])), n).items())

    py.sign_in("al_indigo", "dca63z15bu")
    trace0 = Bar(
        name=str(n)+u'-граммы ридов Класса1 ("хорошие риды")',
        x=[k for (k, v) in oases_distr],
        y=[v for (k, v) in oases_distr],
        opacity=0.9,
        marker=Marker(color='black')

    )

    trace1 = Bar(
        name=str(n)+u'-граммы ридов Класса2 ("плохие риды")',
        x=[k for (k, v) in oases_bad_distr],
        y=[v for (k, v) in oases_bad_distr],
        opacity=0.9,
        marker=Marker(color='grey')
    )

    trace2 = Bar(
        name=str(n)+u'-граммы ридов Класса1 ("хорошие риды")',
        x=[k for (k, v) in trinity_distr],
        y=[v for (k, v) in trinity_distr],
        opacity=0.9,
        marker=Marker(color='black')
    )

    trace3 = Bar(
        name=str(n)+u'-граммы ридов Класса2 ("плохие риды")',
        x=[k for (k, v) in trinity_bad_distr],
        y=[v for (k, v) in trinity_bad_distr],
        opacity=0.9,
        marker=Marker(color='grey')
    )

    data_oases = Data([trace0, trace1])
    data_trinity = Data([trace2, trace3])

#    print("Oases: " + str(n))
#    print wilcoxon(oases_distr.values(), oases_bad_distr.values())
#    print("Trinity:" + str(n))
#    print wilcoxon(trinity_distr.values(), oases_bad_distr.values())

    if n == 2:
        layout_oases = Layout(
            title=u'Вероятности встречаемости для ' + str(n) + u'-грамм в Классе1 и Классе2 для Oases',
            xaxis={'title': str(n) + u'-граммы'},
            yaxis={'title': u'Вероятность'},
            legend=Legend(x=0, y=1)
        )
        layout_trinity = Layout(
            title=u'Вероятности встречаемости для ' + str(n) + u'-грамм в Классе1 и Классе2 для Trinity',
            xaxis={'title': str(n) + u'-граммы'},
            yaxis={'title': u'Вероятность'},
            legend=Legend(x=0, y=1)
        )
    else:
        layout_oases = Layout(
            title=u'Вероятности встречаемости для ' + str(n) + u'-грамм в Классе1 и Классе2 для Oases',
            xaxis={'title': str(n) + u'-граммы', 'autotick': False, 'dtick': int(pow(5, n - 2))},
            yaxis={'title': u'Вероятность'},
            legend=Legend(x=0, y=1)
        )
        layout_trinity = Layout(
            title=u'Вероятности встречаемости для ' + str(n) + u'-грамм в Классе1 и Классе2 для Trinity',
            xaxis={'title': str(n) + u'-граммы', 'autotick': False, 'dtick': int(pow(5, n - 2))},
            yaxis={'title': u'Вероятность'},
            legend=Legend(x=0, y=1)
        )

    fig_oases = Figure(data=data_oases, layout=layout_oases)
    fig_trinity = Figure(data=data_trinity, layout=layout_trinity)

    plot_url1 = py.plot(fig_oases, filename=str(n)+'-grams-histogram-oases')
    time.sleep(10)
    plot_url2 = py.plot(fig_trinity, filename=str(n)+'-grams-histogram-trinity')



#(t, pv) = wilcoxon([float(g)/reference_len for g in reference_distr], [float(g)/all_r_len for g in all_reads_dist])
#print t, pv
