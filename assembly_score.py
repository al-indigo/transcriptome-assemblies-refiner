from classifier import *
from tr_parser import *
from metrics import *


def get_score_for_everything():
    (oases_alignment_data, trinity_alignment_data) = get_alignment_data("data/results_Oases.txt",
                                                                        "data/results_Trinity.txt")

    (ref, oases_reads, oases_name_index, trinity_reads, trinity_name_index) = get_assemblies("data/ref_for_reads.fasta",
                                                                                             "data/Oases.fasta",
                                                                                             "data/Trinity.fasta")

    pdb.set_trace()