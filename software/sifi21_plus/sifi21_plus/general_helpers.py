#!/usr/bin/python

import re
import json
from collections import Counter

def iterparse(j):
    """Work around because json can not load multiple objects.
       Taken from http://stackoverflow.com/questions/22112439/valueerror-extra-data-while-loading-json."""
    nonspace = re.compile(r'\S')
    decoder = json.JSONDecoder()
    pos = 0
    while True:
        matched = nonspace.search(j, pos)
        if not matched:
            break
        pos = matched.start()
        decoded, pos = decoder.raw_decode(j, pos)
        yield decoded


def prepare_json_data(f_in):
    """Prepares the json file for easier parsing."""
    sifi_data = open(f_in, "r").read()
    data = list(iterparse(sifi_data))[0]
    return data

def get_table_data(f_in):
    """Extracts a summary of all and efficient hits."""
    query = prepare_json_data(f_in)
    # Get number of hits per query
    hit_counter = Counter(player['hit_name'] for player in query)
    efficicent_counter = Counter(player['hit_name'] for player in query if player['is_efficient'])
    #hit_overview = hit_counter.most_common()
    table_data = []
    for x in hit_counter.most_common():
        for y in efficicent_counter.most_common():
            if x[0] == y[0]:
                #print y[1]
                table_data.append([x[0], x[1], y[1]])
        if x[0] not in list(efficicent_counter):
            table_data.append([x[0], x[1], 0])
    return table_data
