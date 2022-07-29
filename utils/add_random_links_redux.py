import os
import ipdb
import argparse
import numpy as np
import pandas as pd
import numpy.random as random 

from progressbar import ProgressBar


def read_input(filename):
    """ return list of timestamps, and 2 sets of nodes (can overlap when 
    input graph is not bipartite
    """
    linkstream = 0
    timestamps = set()
    node_set1 = set()
    node_set2 = set()
    with open(filename, 'r') as fin:
        stream = fin.readlines()
        for line in stream:
            t, u, v = line.strip().split()
            #linkstream.append((int(float(t)), u, v))
            linkstream += 1
            timestamps.add(int(float(t)))
            node_set1.add(u)
            node_set2.add(v)
    node_set1 = list(node_set1)
    node_set2 = list(node_set2)
    timestamps = list(timestamps)
    return linkstream, timestamps, node_set1, node_set2

def add_rdm_links(linkstream, timestamps, node_set1, node_set2, p_rdm, output_folder):
    """ add N_rdm random link to link stream"""
    #def check_loop(links):
    #    loop_flag = False
    #    for _, n1, n2 in rdm_links:
    #        if n1 == n2:
    #            loop_flag = True 
    #            break
    #    if loop_flag:
    #        np.random.shuffle(links)
    #        check_loop(links)

    N_rdm = int(p_rdm * 0.01 * linkstream)
    #N_rdm2 = N_rdm
    #if len(node_set1) < 10000:
    #    N_rdm = 2 * N_rdm
    #rdm_nodes1 = random.choice(list(node_set1), size=N_rdm, replace=True)
    #rdm_nodes2 = random.choice(list(node_set2), size=N_rdm, replace=True)
    rdm_nodes1 = []
    rdm_nodes2 = []
    def pick_nodes(node_set1, node_set2):
        same = True
        while same:
            n1 = random.randint(0, len(node_set1))
            n2 = random.randint(0, len(node_set2))
            if node_set1[n1] == node_set2[n2]:
                same = True
            else:
                same = False
        #n1, n2 = random.choice(node_set1, size=2, replace=True)
        
        return node_set1[n1], node_set2[n2]

    #rdm_times = random.choice(list(timestamps), size=N_rdm, replace=True)
    pbar = ProgressBar()
    print(f'creating {N_rdm} random links')
    with open(os.path.join(output_folder, "added_links.txt"), 'w') as fout:#,open(os.path.join(output_folder, "is_fraud.txt"), 'w') as anout : 

        for j in pbar(range(N_rdm)):
            #if j % 100 == 0:
            #    print(f'created {j} random links')
            n1, n2 = pick_nodes(node_set1, node_set2)
            #n1, n2 = random.choice(node_set1, size=2, replace=True)
            #rdm_nodes1.append(n1)
            #rdm_nodes2.append(n2)
            rdm_time = random.randint(0, len(timestamps))
            
            fout.write(f'{timestamps[rdm_time]} {n1} {n2}\n')
            #anout.write(f"{idx},{fraud}\n")

    #timestamps |= set(rdm_times)
    #linkstream += list(zip(rdm_times, rdm_nodes1, rdm_nodes2))
    #rdm_links = list(zip(rdm_times, rdm_nodes1, rdm_nodes2))
    
    # check there's no self loop
    #if (len(node_set1) < 10000):
    #    new_links = []
    #    for t, n1, n2 in rdm_links:
    #        if n1 != n2:
    #            new_links.append((t,n1,n2))
    #    rdm_links = new_links[:N_rdm2]
    #else:
    #    check_loop(rdm_links)
    #return rdm_links

def sort_links(linkstream, rdm_links):
    is_false = [0] * len(linkstream) + [1] * len(rdm_links)
    all_links = list(zip(linkstream + rdm_links, is_false))
    sorted_links = sorted(all_links, key=lambda x: x[0][0])
    return sorted_links

def write_new_link_stream(linkstream, is_false, output_folder):
    with open(os.path.join(output_folder, "added_links.txt"), 'w') as fout,open(os.path.join(output_folder, "is_fraud.txt"), 'w') as anout : 
        anout.write("interaction_id,is_fraud\n")
        for idx, stream in enumerate(linkstream):
            (t, u, v) = stream
            fraud = is_false[idx]
            fout.write(f'{t} {u} {v}\n')
            anout.write(f"{idx},{fraud}\n")

def main():
    parser = argparse.ArgumentParser(description='link stream metrics')
    parser.add_argument('dataset', type=str,
            help='path to the dataset')
    parser.add_argument('Nlinks', type=int,
            help='number of random links to add')
    parser.add_argument('output_folder', type=str,
            help='path to write the output')

    args = parser.parse_args()
    print("reading input")
    ls, t, nodes1, nodes2 = read_input(args.dataset)
    print("adding random links")
    add_rdm_links(ls, t, nodes1, nodes2, args.Nlinks, args.output_folder)
    #is_false = [0] * len(ls) + [1] * len(rdm_links)
    
    #sorted_links = sort_links(ls, rdm_links)
    #write_new_link_stream(sorted_links, args.output_folder)
    #print("writing output")
    #write_new_link_stream(rdm_links, is_false, args.output_folder)


if __name__ == "__main__":
    main()
