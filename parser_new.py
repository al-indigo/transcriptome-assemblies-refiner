from Levenshtein import *

from warnings import warn

import swalign

import pickle

import numpy as np

import pdb

import sys

from sklearn import svm
from sklearn.feature_extraction.text import CountVectorizer


class SymbolTokenizer(object):
    def __call__(self, text):
        tokens = text
        return tokens


def extract_features(vectorizer, texts):
    features = []
    for s in texts:
        features.append(vectorizer.transform(s))
    return features


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
            reads_seq.append(li[0].rstrip('\n'))
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

    top_bound = 0.99
    bottom_bound = 0.30

    for aln in o_res:
        if i > len(o_name_index):
            break

            # print "stage"

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
                # pdb.set_trace()
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
                    # print "success!"
                    break

        i += 1
    return class_trinity_t



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


def train_classifier(class_n, reads_seq, classification_kernel, n_gram_low, n_gram_high, ):
    classifier_t = svm.SVC(kernel=classification_kernel)
    vectorizer = CountVectorizer(ngram_range=(n_gram_low, n_gram_high), tokenizer=SymbolTokenizer())

    corpus_t_1 = class_n
    corpus_t_2 = list(set([item for item in reads_seq if item not in class_n]))
    if len(corpus_t_2) > len(corpus_t_1):
        corpus_t_2 = corpus_t_2[0:len(corpus_t_1)]
    corpus_t = corpus_t_1 + corpus_t_2

    print("Corpus prepared")

    x_t = vectorizer.fit_transform(corpus_t)

    print("Corpus vectorized")

    y_t_1 = [0] * len(corpus_t_1)
    y_t_2 = [1] * len(corpus_t_2)

#    y = y_t + y_o
    y_t = y_t_1 + y_t_2

    classifier_t.fit(x_t, y_t)


    return vectorizer, classifier_t

def classify(vectorizer, classifier, window_size, filename_out, read_name):

    k = 10000
    class_count = 0
    class_uncount = 0
    results_dict = dict()
    for window in range(0, len(reads_seq), k):
        t = reads_seq[window_size:window_size + k]

        #        print("Constructing data set")
        data_set = vectorizer.transform(t)

        #        print("Done; reads count total:", str(len(reads_seq)))
        try:
            classified = classifier.predict(data_set)
        except:
            pdb.set_trace()

            #        print("Predictions ready")

        j = 0
        for i in classified:
            if i == 0:
                class_count += 1
            if i == 1:
                class_uncount += 1
            results_dict[read_name[j]] = i
            j += 1

        f = file(filename_out, "w")
        pickle.dump(results_dict, f)
        f.close()
            #                if (class_t_count + class_o_count) % 10000 == 0:
            #                    print (str(class_t_count + class_o_count), " of ", str(len(reads_seq)))
            #                    print (class_t_count, class_o_count)

    print ("Everything done for ", str(ng), "in", str(int(time.time() - start)), "seconds")
    print (class_count, class_uncount)

    return results_dict

if __name__ == '__main__':
    (reads_n, reads_seq) = get_reads()
    # if False:
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

    class_trinity_t = reads_for_class(o_res, o_name_index, o_dict_t_r, t_dict_r_t, t_dict, t_pairs, t_res, reads_dict)
    f = file("/home/alex/another/diploma/parser/Class_Trinity.txt", "w")
    pickle.dump(class_trinity_t, f)
    f.close()

    class_trinity_o = reads_for_class(t_res, t_name_index, t_dict_t_r, o_dict_r_t, o_dict, o_pairs, o_res, reads_dict)
    f = file("/home/alex/another/diploma/parser/Class_Oases.txt", "w")
    pickle.dump(class_trinity_o, f)
    f.close()

    #    with open("/home/alex/another/diploma/parser/Class_Trinity.txt", 'r') as f:
    #        class_trinity_t = pickle.load(f)

    #    with open("/home/alex/another/diploma/parser/Class_Oases.txt", 'r') as f:
    #        class_trinity_o = pickle.load(f)
    if False:
        intersection = set(class_trinity_t) & set(class_trinity_o)
        print("Classes intersection: ", str(len(intersection)))
        class_trinity_t = [item for item in class_trinity_t if item not in intersection]
        print ("class trinity_t length: ", str(len(class_trinity_t)))
        class_trinity_o = [item for item in class_trinity_o if item not in intersection]
        print("class trinity_o length: ", str(len(class_trinity_o)))
        print("And now intersection is: ", str(len(list(set(class_trinity_t) & set(class_trinity_o)))))

    print reads_seq[1]

    import time

    '''

    for ng in [1, 2, 3, 4, 5, 6, 8, 9, 12, 24]:
        start = time.time()
        classifier = svm.SVC(kernel='rbf')
        vectorizer = CountVectorizer(ngram_range=(ng, ng), tokenizer=SymbolTokenizer())

#        print ("Preparing corpus")

    #    corpus = class_trinity_t + class_trinity_o
        corpus = class_trinity_o + class_trinity_t

#        print("Corpus prepared")

        x = vectorizer.fit_transform(corpus)

#        print("Corpus vectorized")

        y_t = [0] * len(class_trinity_t)
        y_o = [1] * len(class_trinity_o)

    #    y = y_t + y_o
        y = y_o + y_t

        classifier.fit(x, y)

#        print("Training finished")

        k = 10000
        class_t_count = 0
        class_o_count = 0
        for window in range(0, len(reads_seq) - 3200000, k):
            t = reads_seq[window:window+k]

    #        print("Constructing data set")
            data_set = vectorizer.transform(t)

    #        print("Done; reads count total:", str(len(reads_seq)))

            classified = classifier.predict(data_set)

    #        print("Predictions ready")

            for i in classified:
                if i == 0:
                    class_t_count += 1
                if i == 1:
                    class_o_count += 1
#                if (class_t_count + class_o_count) % 10000 == 0:
#                    print (str(class_t_count + class_o_count), " of ", str(len(reads_seq)))
#                    print (class_t_count, class_o_count)

        print ("Everything done for ", str(ng), "in", str(int(time.time() - start)), "seconds")
        print (class_t_count, class_o_count)

    '''

    '''
    for ng in [1, 4]:
        start = time.time()
        classifier_t = svm.SVC(kernel='poly')
#        classifier_o = svm.SVC(kernel='poly')
#        vectorizer = CountVectorizer(ngram_range=(ng, ng), tokenizer=SymbolTokenizer())

        print ("Preparing corpus")
        #    corpus = class_trinity_t + class_trinity_o

    '''

    '''
######################          GRAPHS AND ANOTHER        #####################################
        import plotly.plotly as py
        from plotly.graph_objs import *
        import matplotlib.pyplot as plt

        py.sign_in("aimly", "qpwu7ti8h4")

        def unique_filter(list_to_filter):
            unique_counts = dict()
            for item in list_to_filter:
                if unique_counts.has_key(item):
                    unique_counts[item] += 1
                else:
                    unique_counts[item] = 0
            return unique_counts

        def contain_filter(list_to_filter, where_to_search):
            counts = dict()
            for item in list_to_filter:
                counts[item] = 0
                for read in where_to_search:
                    if item == read:
                        counts[item] += 1
            return counts


        unique_reads = unique_filter(reads_seq)
        unique_trinity = contain_filter(class_trinity_t, reads_seq)
        unique_oases = contain_filter(class_trinity_o, reads_seq)

        def assert_len(reads):
            length = len(reads[0])
            for read in reads:
                if len(read) != length:
                    pdb.set_trace()
            print length

        assert_len(reads_seq)
        assert_len(class_trinity_t)
        assert_len(class_trinity_o)

        #        sorted_x = sorted(unique_reads, key=unique_reads.get, reverse=True)

        y_reads = sorted(unique_reads.values(), reverse=True)
        x_reads = range(len(unique_reads))

        y_trinity = sorted(unique_trinity.values(), reverse=True)
        x_trinity = range(len(unique_trinity))

        y_oases = sorted(unique_oases.values(), reverse=True)
        x_oases = range(len(unique_oases))

        trace1 = Bar(x=x_reads, y=y_reads, name="Reads full")
        trace2 = Bar(x=x_trinity, y=y_trinity, name="Reads trinity")
        trace3 = Bar(x=x_oases, y=y_oases, name="Reads oasis")

        width = 0.01

        #rect_reads = plt.bar(x_tr, y_reads, width, color='b', label='Reads full')

        #plt.show()
        #sleep(15)
        #plt.close()
        data = Data([trace2, trace3])
    #        plot_url = py.plot(data, filename='test')
    '''
############################################################################
    '''
        corpus_t_1 = class_trinity_t
        corpus_t_2 = list(set([item for item in reads_seq if item not in class_trinity_t]))
        if len(corpus_t_2) > len(corpus_t_1):
            corpus_t_2 = corpus_t_2[0:len(corpus_t_1)]
        corpus_t = corpus_t_1 + corpus_t_2

        print("Corpus prepared")

        x_t = vectorizer.fit_transform(corpus_t)

        print("Corpus vectorized")

        y_t_1 = [0] * len(corpus_t_1)
        y_t_2 = [1] * len(corpus_t_2)

    #    y = y_t + y_o
        y_t = y_t_1 + y_t_2

        classifier_t.fit(x_t, y_t)
    '''
    '''
    corpus_o_1 = class_trinity_o
    corpus_o_2 = list(set([item for item in reads_seq if item not in class_trinity_o]))
    if len(corpus_o_2) > len(corpus_o_1):
        corpus_o_2 = corpus_o_2[0:len(corpus_o_1)]
    corpus_o = corpus_o_1 + corpus_o_2

    print("Corpus prepared")

    x_o = vectorizer.fit_transform(corpus_o)

    print("Corpus vectorized")

    y_o_1 = [0] * len(corpus_o_1)
    y_o_2 = [1] * len(corpus_o_2)

    #    y = y_t + y_o
    y_o = y_o_1 + y_o_2

    classifier_o.fit(x_o, y_o)

    print("Training finished")

    k = 10000
    class_t_count = 0
    class_t_uncount = 0
    class_o_count = 0
    class_o_uncount = 0
    for window in range(0, len(reads_seq) - 3100000, k):
        t = reads_seq[window:window + k]

        #        print("Constructing data set")
        data_set_t = vectorizer.transform(t)
        data_set_o = vectorizer.transform(t)

        #        print("Done; reads count total:", str(len(reads_seq)))
        try:
            classified_t = classifier_t.predict(data_set)
            classified_o = classifier_o.predict(data_set)
        except:
            pdb.set_trace()

            #        print("Predictions ready")

        for i in classified_t:
            if i == 0:
                class_t_count += 1
            if i == 1:
                class_t_uncount += 1
        for j in classified_o:
            if j == 0:
                class_o_count += 1
            if j == 1:
                class_o_uncount += 1
            #                if (class_t_count + class_o_count) % 10000 == 0:
            #                    print (str(class_t_count + class_o_count), " of ", str(len(reads_seq)))
            #                    print (class_t_count, class_o_count)

    print ("Everything done for ", str(ng), "in", str(int(time.time() - start)), "seconds")
    print (class_t_count, class_t_uncount)
    print (class_o_count, class_o_uncount)
    '''
    (vectorizer_o, classifier_o) = train_classifier(class_trinity_o, reads_seq, "rbf", 1, 4)
    (vectorizer_t, classifier_t) = train_classifier(class_trinity_t, reads_seq, "rbf", 1, 4)
    (results_dict_o) = classify(vectorizer_o, classifier_o, 10000, "/home/alex/another/diploma/parser/oases_classified", reads_n)
    (results_dict_t) = classify(vectorizer_t, classifier_t, 10000, "/home/alex/another/diploma/parser/trinity_classified", reads_n)

    # import cProfile

    #cProfile.run('main()')