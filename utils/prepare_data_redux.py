import os
import ipdb
import argparse
import datetime
import numpy as np
import pandas as pd
import numpy.random as random 

from progressbar import ProgressBar

def map_stream(stream, output_stream, output_mapping):
    """ Assume sorted  stream in input, integer timestamps"""

    node_map = dict()
    node_idx = 0
    pbar = ProgressBar()
    #def map_node(u):
    #    if u in node_map:
    #        u_mapped = node_map[u]
    #    else:
    #        u_mapped = node_idx
    #        node_idx += 1
    #    return u_mapped

    with open(stream, 'r') as fin, open(output_stream, 'w') as fout:
        for line_idx, line in enumerate(fin):

            _t,u,v = line.strip().split(',')
            t = datetime.datetime.fromisoformat(_t).timestamp()
            if line_idx == 0:
                t0 = int(t)
            if u in node_map:
                u_mapped = node_map[u]
            else:
                u_mapped = node_idx
                node_map[u] = u_mapped

                node_idx += 1
            if v in node_map:
                v_mapped = node_map[v]
            else:
                v_mapped = node_idx
                node_map[v] = v_mapped
                node_idx += 1

            #u_mapped = map_node(u)
            #v_mapped = map_node(v)
            fout.write(f'{int(t) - t0} {u_mapped} {v_mapped}\n')
    with open(output_mapping, 'w') as fout:
        for key in node_map:
            fout.write(f'{key} {node_map[key]}\n')

def main():
    parser = argparse.ArgumentParser(description='link stream metrics')
    parser.add_argument('dataset', type=str,
            help='path to the dataset')
    parser.add_argument('output_stream', type=str,
            help='path to output stream')
    parser.add_argument('output_mapping', type=str,
            help='path to output mapping')

    args = parser.parse_args()
    map_stream(args.dataset, args.output_stream, args.output_mapping)

if __name__ == "__main__":
    main()
