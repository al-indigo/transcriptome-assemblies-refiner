import pickle
from sklearn import svm, ensemble, semi_supervised, covariance
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

import pdb
def train_classifier(class_n, reads_seq, method, n_gram_low, n_gram_high):
    if 'svc' in method:
        if 'count' in method:
            classifier_obj = svm.SVC(kernel=method['svc'], class_weight=None, degree=method['count'])
        else:
            classifier_obj = svm.SVC(kernel=method['svc'], class_weight=None)
    elif 'rfc' in method:
        classifier_obj = ensemble.RandomForestClassifier(n_estimators=method[method.keys()[0]], min_samples_split=2, min_samples_leaf=1)
    elif 'oneclasssvm' in method:
        classifier_obj = svm.OneClassSVM(max_iter=100)
    elif 'covariance' in method:
        classifier_obj = covariance.EllipticEnvelope()
    else:
        classifier_obj = semi_supervised.LabelSpreading(kernel='knn')

    vectorizer_obj = CountVectorizer(ngram_range=(n_gram_low, n_gram_high), tokenizer=SymbolTokenizer())

    corpus_good = class_n
    corpus_bad = list(set([item for item in reads_seq if item not in class_n]))
    if len(corpus_bad) > len(corpus_good):
        corpus_unlabeled = corpus_bad[len(corpus_good):len(corpus_good)*2]
        corpus_bad = corpus_bad[0:len(corpus_good)]

    corpus = corpus_good + corpus_bad #+ corpus_unlabeled

#    print("Corpus prepared")
    x = vectorizer_obj.fit_transform(corpus).todense()
    y = [1] * len(corpus_good) + [0] * len(corpus_bad) #+ [-1] * len(corpus_unlabeled)
    #x = vectorizer_obj.fit_transform(corpus_good)
#    print("Corpus vectorized")

    classifier_obj.fit(x, y)
    #classifier_obj.fit(x)

    return vectorizer_obj, classifier_obj


def classify(vectorizer_obj, classifier_obj, window_size, filename_out, reads_names, reads_seq):
    class_good = 0
    class_bad = 0
    f = file(filename_out, "wb")
    for pos in range(0, len(reads_seq), window_size):
        results_dict = dict()
        if pos + window_size < len(reads_seq):
            t = reads_seq[pos:pos + window_size]
        else:
            t = reads_seq[pos:len(reads_seq)]
        data_set = vectorizer_obj.transform(t).todense()
        classified = classifier_obj.predict(data_set)

        j = 0
        for i in classified:
#            if i == +1:
            if i == 0:
                class_good += 1
            if i == 1:
#            if i == -1:
                class_bad += 1
            results_dict[reads_names[pos + j]] = i
            j += 1
        pickle.dump(results_dict, f)

    print (class_good, class_bad, filename_out)
    f.close()

    f = file(filename_out + "-SUMMARY", "w")
    f.write(filename_out + "\t" + str(class_good) + "\t" + str(class_bad))
    f.close()


    return