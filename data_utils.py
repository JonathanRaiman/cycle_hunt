import gzip
import numpy as np
import dali.core as D

from dali.data.batch import create_lines_batch
from xml_cleaner     import tokenize

def to_words(txt):
    return " ".join([word for sentence in tokenize(txt) for word in sentence])

class Example(object):
    def __init__(self, text, triggers):
        self.text = to_words(text)
        self.triggers = triggers

def load_dataset(path):
    examples = {}
    with gzip.open(path, "rt") as f:
        for l in f:
            example, *labels = l.strip().split("\t")
            ex = Example(example, labels)
            for label in labels:
                if label in examples:
                    examples[label].append(ex)
                else:
                    examples[label] = [ex]
    return examples

def produce_negative_example_set(positives, labels_vocab, negatives):
    batch_size = len(positives) + negatives
    keys       = np.zeros(batch_size, dtype=np.int32)
    labels     = np.zeros(batch_size, dtype=np.float32)

    keys[:len(positives)] = positives
    labels[:len(positives)] = 1.0

    # sample some element from the vocab:
    if negatives > 0:
        samples = np.random.choice(len(labels_vocab), size=batch_size, replace=False)
        positives_set = set(positives)
        idx = len(positives)
        for sample in samples:
            if sample not in positives_set:
                keys[idx] = sample
            idx += 1
            if idx == batch_size:
                break
    return D.Mat(keys, borrow=True, dtype=np.int32), D.Mat(labels, borrow=True, dtype=np.float32)

def example_to_data(example, word_vocab, labels_vocab, batch_size):
    x = create_lines_batch([example.text], word_vocab)
    keys, labels = produce_negative_example_set(
        labels_vocab.encode(
            [trig for trig in example.triggers if trig in labels_vocab]
        ),
        labels_vocab,
        batch_size
    )
    return x, keys, labels
