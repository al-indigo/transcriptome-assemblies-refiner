from classifier import *
from tr_parser import *
from metrics import *

import multiprocessing

import pdb
import operator


def classification_cycle(kernel,
                         ngram,
                         class_good_oases,
                         class_good_trinity,
                         reads_names,
                         reads_seq,
                         fileout_suffix,
                         is_threaded):
    #print ("thread started")
    (vectorizer_o, classifier_o) = train_classifier(class_good_oases, reads_seq, kernel, ngram, ngram)
    (vectorizer_t, classifier_t) = train_classifier(class_good_trinity, reads_seq, kernel, ngram, ngram)

    classify(vectorizer_obj=vectorizer_o,
             classifier_obj=classifier_o,
             window_size=10000,
             filename_out="results/oases_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix,
             reads_names=reads_names,
             reads_seq=reads_seq)

    classify(vectorizer_obj=vectorizer_t,
             classifier_obj=classifier_t,
             window_size=10000,
             filename_out="results/trinity_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix,
             reads_names=reads_names,
             reads_seq=reads_seq)


def validation_cycle(kernel,
                     ngram,
                     class_good_oases,
                     class_good_trinity,
                     reads_names,
                     reads_seq,
                     fileout_suffix,
                     is_threaded):
    #print ("thread started")
    if 'count' in kernel:
        fileout_suffix += '-degree-' + str(kernel['count']) + '-'
    oases_training_set = class_good_oases[0:int(len(class_good_oases)*0.75)]
    trinity_training_set = class_good_trinity[0:int(len(class_good_trinity)*0.75)]
    oases_test_set = class_good_oases[int(len(class_good_oases)*0.75)+1:int(len(class_good_oases)-1)]
    trinity_test_set = class_good_trinity[int(len(class_good_trinity)*0.75)+1:int(len(class_good_trinity)-1)]

    (vectorizer_o, classifier_o) = train_classifier(oases_training_set, reads_seq, kernel, 4, ngram)
    (vectorizer_t, classifier_t) = train_classifier(trinity_training_set, reads_seq, kernel, 4, ngram)

    if is_threaded:
        process_list = []
        p = multiprocessing.Process(target=classify, args=(vectorizer_o,
                                                           classifier_o,
                                                           10000,
                                                           "validation/oases_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-fit",
                                                           reads_names,
                                                           oases_test_set))
        process_list.append(p)
        p.start()

        p = multiprocessing.Process(target=classify, args=(vectorizer_t,
                                                           classifier_t,
                                                           10000,
                                                           "validation/trinity_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-fit",
                                                           reads_names,
                                                           trinity_test_set))
        process_list.append(p)
        p.start()

        test_list = list(operator.xor(set(reads_seq), set(class_good_oases)))

        p = multiprocessing.Process(target=classify, args=(vectorizer_o,
                                                           classifier_o,
                                                           10000,
                                                           "validation/oases_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-not-fit",
                                                           reads_names,
                                                           test_list[0:10000]))
        process_list.append(p)
        p.start()

        test_list = list(operator.xor(set(reads_seq), set(class_good_trinity)))
        p = multiprocessing.Process(target=classify, args=(vectorizer_t,
                                                           classifier_t,
                                                           10000,
                                                           "validation/trinity_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-not-fit",
                                                           reads_names,
                                                           test_list[0:10000]))
        process_list.append(p)
        p.start()

        for p in process_list:
            p.join()
    else:
        classify(vectorizer_o,
                 classifier_o,
                 10000,
                 "validation/oases_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(
                     kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-fit",
                 reads_names,
                 oases_test_set)

        classify(vectorizer_t,
                 classifier_t,
                 10000,
                 "validation/trinity_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-fit",
                 reads_names,
                 trinity_test_set)

        test_list = list(operator.xor(set(reads_seq), set(class_good_oases)))
        classify(vectorizer_o,
                 classifier_o,
                 10000,
                 "validation/oases_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-not-fit",
                 reads_names,
                 test_list[0:10000])

        test_list = list(operator.xor(set(reads_seq), set(class_good_trinity)))
        classify(vectorizer_t,
                 classifier_t,
                 10000,
                 "validation/trinity_classified-" + str(ngram) + "-" + kernel.keys()[0] + "-" + str(kernel[kernel.keys()[0]]) + "-" + fileout_suffix + "-must-not-fit",
                 reads_names,
                 test_list[0:10000])


def launch_classification_threads(class_good_oases,
                                  class_good_trinity,
                                  reads_names,
                                  reads_seq,
                                  purpose,
                                  fileout_suffix,
                                  is_threaded=False):
    target = None

    if purpose == "real-life":
        target = classification_cycle
    elif purpose == "validate":
        print ("cross-validating")
        fileout_suffix += "-validating"
        target = validation_cycle

#    for ngram in [2, 3, 4, 5, 6, 7, 8]:
#    kernels = [{'svc': 'rbf'}, {'svc': 'poly'}, {'svc': 'sigmoid'}, {'rfc': 2}, {'rfc': 4}, {'rfc': 10}, {'rfc': 100}, {'rfc': 200}, {'rfc': 1000}]
    kernels = [{'rfc': 1000}]
#    kernels = [{'svc': 'poly', 'count': 3},
#               {'svc': 'poly', 'count': 4},
#               {'svc': 'poly', 'count': 5},
#               {'svc': 'poly', 'count': 6},
#               {'svc': 'poly', 'count': 7},
#               {'svc': 'poly', 'count': 8}]
#    kernels = [{'lp': ''}]
#    kernels = [{'oneclasssvm': ''}]
#    kernels = [{'covariance': ''}]

    for kernel in kernels:
        for ngram in [6]:
            target(kernel,
                   ngram,
                   class_good_oases,
                   class_good_trinity,
                   reads_names,
                   reads_seq,
                   fileout_suffix,
                   is_threaded)



def main():
    compute_for_real = True
    validate = False
    for bottom_bound in [0.3]:
        #for top_bound in [0.99, 0.95, 0.9]:
        for top_bound in [0.9]:
            fileout_suffix = "-" + str(bottom_bound) + "-" + str(top_bound) + "-"

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

            if compute_for_real:
                launch_classification_threads(class_good_oases,
                                              class_good_trinity,
                                              reads_names,
                                              reads_seq,
                                              "real-life",
                                              fileout_suffix)
            if validate:
                launch_classification_threads(class_good_oases,
                                              class_good_trinity,
                                              reads_names,
                                              reads_seq,
                                              "validate",
                                              fileout_suffix)


start = time.time()

main()

print ('Completed in ' + str(int(time.time() - start)) + " seconds")


    # import cProfile

    #cProfile.run('main()')