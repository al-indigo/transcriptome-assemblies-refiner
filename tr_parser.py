import pickle


###################     READS FROM FILE       ########################
def get_reads(filename):
    f = file(filename, "r")

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


######################      ALIGNMENTS    #########################
def get_alignment_data(oases_alignment_filename, trinity_alignment_filename):
    oases_alignment_data = []
    with open(oases_alignment_filename, "r") as f:
        oases_alignment_data = pickle.load(f)
        f.close()

    trinity_alignment_data = []
    with open(trinity_alignment_filename, "r") as f:
        trinity_alignment_data = pickle.load(f)
        f.close()

    return oases_alignment_data, trinity_alignment_data


#################    READS FROM SAM FILE #####################
def get_reads_for_assembler(sam_reads_filename):
    f = file(sam_reads_filename, "r")

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


#########################   ASSEMBLIES READING   ###########################
#TODO: implement it correctly (this implementation is dumb)
def get_assemblies(reference_filename, oases_filename, trinity_filename):
    f = open(reference_filename, "r")
    ref_transcripts = []

    tmp = ""
    while 1:
        l = f.readline()
        if not l:
            ref_transcripts.append(tmp)
            tmp = ""
            break
        if l[0] != '>':
            tmp += l[:-1]

        else:
            ref_transcripts.append(tmp)
            tmp = ""

    f.close()

    f = open(oases_filename, "r")
    oases_transcripts = []
    oases_name_index = []

    tmp = ""
    while 1:
        l = f.readline()
        if not l:
            oases_transcripts.append(tmp)
            tmp = ""
            break
        if l[0] != '>':
            tmp += l[:-1]
        else:
            oases_transcripts.append(tmp)
            tmp = ""
            massive = l.split(' ')
            oases_name_index.append(massive[0][1:-1])

    f.close()

    f = open(trinity_filename, "r")
    trinity_transcripts = []
    trinity_name_index = []

    tmp = ""
    while 1:
        l = f.readline()
        if not l:
            trinity_transcripts.append(tmp)
            tmp = ""
            break
        if l[0] != '>':
            tmp += l[:-1]
        else:
            trinity_transcripts.append(tmp)
            tmp = ""
            massive = l.split(' ')
            trinity_name_index.append(massive[0][1:])

    f.close()
    return ref_transcripts, oases_transcripts, oases_name_index, trinity_transcripts, trinity_name_index


def make_index_reads_to_transcripts(reads_names, transcripts_names):
    index = dict()

    j = 0
    for i in reads_names:
        index[i] = transcripts_names[j]
        j += 1

    return index


def make_index_transcripts_to_reads(name_index, transcripts_names, reads_names):
    index = dict()

    for i in name_index:
        index[i] = []

    j = 0
    t = 0

    for i in transcripts_names:
        if index.has_key(i):
            index[i].append(reads_names[j])
        j += 1

    return index


def make_index_by_name(name_index):
    index = dict()

    j = 0
    for i in name_index:
        index[i] = j
        j += 1

    return index