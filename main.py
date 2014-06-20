from classifier import *
from tr_parser import *
from metrics import *

from multiprocessing import Process

import pdb


def classification_cycle(kernel,
                         ngram,
                         class_good_oases,
                         class_good_trinity,
                         reads_names,
                         reads_seq):
    print ("thread started")
    (vectorizer_o, classifier_o) = train_classifier(class_good_oases, reads_seq, kernel, ngram, ngram)
    (vectorizer_t, classifier_t) = train_classifier(class_good_trinity, reads_seq, kernel, ngram, ngram)
    classify(vectorizer_o, classifier_o, 10000, "results/oases_classified-" + str(ngram) + "-" + kernel, reads_names, reads_seq)
    classify(vectorizer_t, classifier_t, 10000, "results/trinity_classified-" + str(ngram) + "-" + kernel, reads_names, reads_seq)


def main():
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
                                         reference_reads_to_transcripts)
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
                                       reference_reads_to_transcripts)
    f = file("data/Class_Oases.txt", "w")
    pickle.dump(class_good_oases, f)
    f.close()

    process_list = []

    for kernel in ['poly', 'rbf']:
        for ngram in [3, 4]:
            p = Process(target=classification_cycle, args=(kernel,
                                                           ngram,
                                                           class_good_oases,
                                                           class_good_trinity,
                                                           reads_names,
                                                           reads_seq))
            process_list.append(p)
            p.start()

    for p in process_list:
        p.join()


main()
    # import cProfile

    #cProfile.run('main()')