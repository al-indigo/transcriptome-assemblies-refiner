import csv, pdb
import pprint
def check_validation_results(assembler):
    missing = 0
    existing = 0
    good = 0
    good_list = []
    kernels = [{'svc': 'rbf'}, {'svc': 'poly'}, {'svc': 'sigmoid'}, {'rfc': 2}, {'rfc': 4}, {'rfc': 10}, {'rfc': 100}, {'rfc': 200}, {'rfc': 1000}]
    for bottom_bound in [0.3, 0.4, 0.5, 0.6]:
        for top_bound in [0.99, 0.95, 0.9]:
            for ngram in [2, 3, 4, 5, 6, 7, 8]:
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
#    pp.pprint(good_list)
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

from tr_parser import get_assemblies
(ref, oases_reads, oases_name_index, trinity_reads, trinity_name_index) = get_assemblies("data/ref_for_reads.fasta",
                                                                                         "data/Oases.fasta",
                                                                                         "data/Trinity.fasta")
from ngram import NGram

n = NGram(N=4, pad_len=0)
grams = dict()
for transcript in ref:
    if transcript == '':
        continue
    ngram_list = list(n._split(transcript))
    for ng in ngram_list:
        if ng == 'TTSG':
            pdb.set_trace()
        if ng in grams:
            grams[ng] += 1
        else:
            grams[ng] = 1

pp = pprint.PrettyPrinter()
pp.pprint(grams)
print len(grams)
