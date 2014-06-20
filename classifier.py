import pickle
from sklearn import svm
from sklearn.feature_extraction.text import CountVectorizer

import time


class SymbolTokenizer(object):
    def __call__(self, text):
        tokens = text
        return tokens


def extract_features(vectorizer, texts):
    features = []
    for s in texts:
        features.append(vectorizer.transform(s))
    return features


def train_classifier(class_n, reads_seq, classification_kernel, n_gram_low, n_gram_high):
    classifier_obj = svm.SVC(kernel=classification_kernel)
    vectorizer_obj = CountVectorizer(ngram_range=(n_gram_low, n_gram_high), tokenizer=SymbolTokenizer())

    corpus_good = class_n
    corpus_bad = list(set([item for item in reads_seq if item not in class_n]))
    if len(corpus_bad) > len(corpus_good)*3:
        corpus_bad = corpus_bad[0:len(corpus_good)*3]
    corpus = corpus_good + corpus_bad

    print("Corpus prepared")

    x = vectorizer_obj.fit_transform(corpus)
    y = [0] * len(corpus_good) + [1] * len(corpus_bad)

    print("Corpus vectorized")

    classifier_obj.fit(x, y)

    return vectorizer_obj, classifier_obj


def classify(vectorizer_obj, classifier_obj, window_size, filename_out, read_name, reads_seq):
    class_good = 0
    class_bad = 0
    f = file(filename_out, "wb")
    for pos in range(0, len(reads_seq), window_size):
        results_dict = dict()
        if pos + window_size < len(reads_seq):
            t = reads_seq[pos:pos + window_size]
        else:
            t = reads_seq[pos:len(reads_seq)]
        data_set = vectorizer_obj.transform(t)
        classified = classifier_obj.predict(data_set)

        j = 0
        for i in classified:
            if i == 0:
                class_good += 1
            if i == 1:
                class_bad += 1
            results_dict[read_name[j]] = i
            j += 1
        pickle.dump(results_dict, f)

    print (class_good, class_bad, filename_out)
    f.close()

    return