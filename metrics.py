import pickle
import pdb


######################      DISTANCES     #######################
def get_distances(similar_transkripts_oases_file,
                  similar_transkripts_trinity_file):
    distances_pairs_for_oases = []
    with open(similar_transkripts_oases_file, "r") as f:
        distances_pairs_for_oases = pickle.load(f)
        f.close()

    distances_pairs_for_trinity = []
    with open(similar_transkripts_trinity_file, "r") as f:
        distances_pairs_for_trinity = pickle.load(f)
        f.close()

    return distances_pairs_for_oases, distances_pairs_for_trinity


def reads_for_class(o_res, o_name_index, o_dict_t_r, t_dict_r_t, t_dict, t_pairs, t_res, reads_dict):
    good_reads = []
    i = 0

    top_bound = 0.99
    bottom_bound = 0.30

    for aln in o_res:
        if i > len(o_name_index):
            break

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
                    good_reads.append(reads_dict[read])
                    # print "success!"
                    break

        i += 1
    return good_reads