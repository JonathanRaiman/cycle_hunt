import os
import gzip
import time
import json
import gzip
import sqlite3
import argparse
import xml_cleaner
import wikipedia_ner.parse

import numpy as np

from multiprocessing import Pool, Queue, Process
from collections     import deque
from os.path         import join, isdir, exists

def uber_decode(s):
    try:
        return s.decode("utf-8")
    except:
        return s.decode("unicode-escape")

def preprocess_title(s):
    return " ".join([word
              for sentence in xml_cleaner.tokenize(s)
              for word in sentence])

def save_graph_from_root(root, graph, outpath, connector="->", max_log1p_degree=7, max_depth = 99999):
    save_nodes = set()
    node_queue = deque()
    node_queue.appendleft(root)
    save_nodes.add(root)
    depth = {root: 0}

    def acceptable_node(node):
        if depth[node] > max_depth:
            return False
        return (
            not node.startswith("Category:")
            or node in graph) and not node.lower().startswith("wikipedia:") and \
            not node.lower().startswith("file:") and \
            not node.lower().startswith("template:")

    with gzip.open(outpath, "wt") as fout:
        while len(node_queue) > 0:
            node = node_queue.pop()
            first = True
            if node in graph:
                if np.log1p(len(graph[node])) <= max_log1p_degree:
                    for child in graph[node]:
                        if child not in depth:
                            depth[child] = depth[node] + 1
                        if acceptable_node(child):
                            if first:
                                fout.write("%s%s%s\n" % (node, connector, child))
                            else:
                                fout.write("%s\n" % (child,))
                            first = False
                            if child not in save_nodes:
                                node_queue.appendleft(child)
                                save_nodes.add(child)
    return save_nodes

def load_index2target(path):
    targets = {}
    with gzip.open(path, "rb") as f:
        for k, l in enumerate(f):
            l2 = uber_decode(l)
            targets[l2.strip()] = k
    return targets

def filter_articles(index2target, articles):
    detected      = set()
    missing       = set()
    iffy_detected = set()
    for title in articles:
        if title in index2target:
            detected.add((title, index2target[title]))
        elif preprocess_title(title) in index2target:
            detected.add((preprocess_title(title), index2target[preprocess_title(title)]))
        else:
            if "(" in title:
                if title.split("(", 1)[0].strip() in index2target:
                    iffy_detected.add((title.split("(", 1)[0].strip(), index2target[title.split("(", 1)[0].strip()]))
                else:
                    missing.add(title)
            else:
                missing.add(title)
    return detected, iffy_detected, missing

def acquire_objs(objnum):
    objs = get_lines_from_db(objnum)
    if objs is not None and type(objs[0]) is not list:
        return objs[0]
    return None

def acquire_objs_worker(jobqueue, outqueue, db_path):
    sqlite_conn = sqlite3.connect(db_path, detect_types=sqlite3.PARSE_DECLTYPES)
    try:
        insert_into_db, update_in_db, update_lines_in_db, get_obj_from_db, get_lines_from_db = wikipedia_ner.parse.sqlite_utils.create_schema(
            sqlite_conn, [
                ("lines", "pickle"),
                ("parents", "pickle")
            ], "articles"
        )
        while True:
            job = jobqueue.get()
            if job is None:
                break
            objs = get_lines_from_db(job)
            if objs is not None and type(objs[0]) is not list:
                corpus = objs[0]
                examples = []
                for example in corpus.example:
                    examples.append(
                        (
                            " ".join(example.words),
                            [trigger.trigger for trigger in example.trigger]
                        )
                    )
                outqueue.put(examples)
            else:
                outqueue.put(None)
    finally:
        sqlite_conn.close()

def collect_corpus(db_path, articles, total=100000, processes=8):
    corpuses = []
    missing = 0
    jobqueue = Queue()
    outqueue = Queue()
    workers = [Process(target=acquire_objs_worker, daemon=True, args=(jobqueue, outqueue, db_path)) for i in range(processes)]
    for worker in workers:
        worker.start()
    examples = []

    submittable = min(len(articles), total)
    submitted = 0
    received  = 0

    def flush_outqueue(received, submitted, missing, batch):
        start = received
        while received < submitted:
            out = outqueue.get()
            received += 1
            if out is not None:
                examples.extend(out)
            else:
                missing += 1
            if received % 100 == 0:
                print("%.3f%% done\r" % (100.0 * received / submittable,), flush=True, end="")
            if received - start >= batch:
                break
        return received, missing

    try:
        for title, num in list(articles)[:total]:
            jobqueue.put(num)
            submitted += 1
            if submitted % 100 == 0:
                #print("%.3f%% submitted\r" % (100.0 * (submitted / ),), flush=True, end="")
                received, missing = flush_outqueue(received, submitted, missing, 100)

        print(" " * 100)

        for worker in workers:
            jobqueue.put(None)

        while received < submitted:
            out = outqueue.get()
            received += 1
            if out is not None:
                examples.extend(out)
            else:
                missing += 1
            if received % 100 == 0:
                print("%.3f%% done\r" % (100.0 * received / submittable,), flush=True, end="")
        print(" " * 100)
    finally:
        for worker in workers:
            worker.join()
    return examples, missing

def vprint(statement, verbosity):
    if verbosity > 0: print(statement)

def filter_triggers(triggers, whitelist):
    kept_triggers = []
    for trigger in triggers:
        if trigger.lower().startswith("category:") or trigger.lower().startswith("category :"):
            continue
        elif len(trigger.strip()) > 0:
            if trigger.strip() in whitelist:
                kept_triggers.append(trigger)
            else:
                new_trigger = (trigger
                    .replace("( ", "(")     \
                    .replace(" )", ")")     \
                    .replace(" , ", ", ")   \
                    .replace(" : ", ":")    \
                    .replace(" 's ", "'s")  \
                    .replace("s ' ", "s' ") \
                    .replace(" - ", "-")    \
                    .replace(" l' ", " l'") \
                    .replace("l ' ", "l' ") \
                    .replace(" d' ", " d'"))
                new_trigger = new_trigger[0].upper() + new_trigger[1:]
                if new_trigger in whitelist:
                    kept_triggers.append(new_trigger)
    return kept_triggers

def save_corpus(examples, whitelist, path, max_labels=None):
    if max_labels is None:
        max_labels = 999999
    with gzip.open(path, "wt") as fout:
        for k, (example, triggers) in enumerate(examples):
            kept_triggers = filter_triggers(triggers, whitelist)
            if len(kept_triggers) > 0 and len(kept_triggers) <= max_labels:
                fout.write(example)
                for trig in kept_triggers:
                    if "Category:" in trig:
                        pass
                    fout.write("\t")
                    fout.write(trig)
                fout.write("\n")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", "--wiki_db",
                        type=str,
                        default="/Users/jonathanraiman/Desktop/datasets/enwiki.db")
    parser.add_argument("-g", "--graph",
                        type=str,
                        default="/Users/jonathanraiman/Desktop/datasets/category_graph.json")
    parser.add_argument("-o", "--output",
                        type=str,
                        required=True)
    parser.add_argument("-i", "--index2target",
                        type=str,
                        default="/Users/jonathanraiman/Desktop/datasets/index2target.clean.gz")
    parser.add_argument("-r", "--root",
                        type=str,
                        required=True)
    parser.add_argument("-v", "--verbosity",
                        type=int, default=1)
    parser.add_argument("-m", "--max_labels", type=int, default=9999999)
    parser.add_argument("--max_log1p_degree", type=int, default=7)
    parser.add_argument("--max_depth",        type=int, default=7)
    parser.add_argument("-t", "--total",      type=int, default=9999999)
    parser.add_argument("-j", "--processes",  type=int, default=8)
    return parser.parse_args()

def load_graph(path):
    with open(path, "rt") as f:
        graph = json.load(f)
    return graph

def validate_arguments(args):
    assert(args.root.startswith("Category:")), \
        "Root must be a category (start with `Category:`)"
    assert((isdir(args.output) or not exists(args.output))), "Output must be a directory."
    assert(args.processes > 0), "Number of processes must be strictly positive."
    os.makedirs(args.output, exist_ok=True)

def main():
    t0 = time.time()
    args = parse_args()
    validate_arguments(args)

    vprint("Loading graph.", args.verbosity)
    graph        = load_graph(args.graph)
    vprint("Done.", args.verbosity)

    vprint("Constructing graph starting at \"" + args.root + "\"", args.verbosity)
    tcollectgraphbegin = time.time()
    visited_nodes = save_graph_from_root(args.root,
                         graph=graph,
                         outpath=join(args.output, "ontology.txt.gz"),
                         connector="->",
                         max_log1p_degree=args.max_log1p_degree,
                         max_depth=args.max_depth
    )
    tcollectgraph = time.time() - tcollectgraphbegin
    print("Collect graph time %.3fs" % (tcollectgraph,))
    vprint("Done.", args.verbosity)

    vprint("Filtering articles using index2target list.", args.verbosity)
    tfilterbegin = time.time()
    detected_articles, iffy_detected_articles, missing = filter_articles(
        index2target=load_index2target(args.index2target),
        articles=(vis for vis in visited_nodes if not vis.startswith("Category:"))
    )
    articles_present = detected_articles | iffy_detected_articles
    vprint("Done. %d missing, %d found (%d exact, %d approximate)" % (
        len(missing), len(articles_present), len(detected_articles), len(iffy_detected_articles)
    ), args.verbosity)
    tfilter = time.time() - tfilterbegin
    print("Filter articles time %.3fs" % (tfilter,))

    vprint("Collecting corpus of examples from db \"%s\"" % (args.wiki_db), args.verbosity)
    tbegincorpus = time.time()
    examples, missing = collect_corpus(
        db_path=args.wiki_db,
        articles=articles_present,
        total=args.total,
        processes=args.processes
    )
    vprint("Done. Got %d examples (%d missing)" % (len(examples), missing), args.verbosity)
    tcollectcorpus = time.time() - tbegincorpus
    print("Collect corpus time %.3fs" % (tcollectcorpus,))

    train_path = join(args.output, "train.tsv.gz")
    vprint("Saving training data under \"%s\"" % (train_path,), args.verbosity)
    save_corpus(
        examples=examples,
        whitelist=visited_nodes,
        path=train_path,
        max_labels=args.max_labels
    )
    vprint("Done", args.verbosity)
    tend = time.time()
    print("Total time %.3fs" % (tend - t0,))

if __name__ == "__main__":
    main()
