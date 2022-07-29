import os
import sys
import ipdb
import argparse
import numpy as np

from contextlib import ExitStack
from sklearn.ensemble import RandomForestClassifier #, AdaBoostClassifier
from sklearn.metrics import f1_score, recall_score, accuracy_score, precision_score, roc_auc_score, mean_squared_error, roc_curve, plot_roc_curve

columns0 = [
        #"number of links",
        #"number of nodes",
        #"total weight",
        #"count nodes degree 1",
        #"count nodes degree 2",
        "degree u",
        "degree v",

        #"max degree",
        #"max weighted degree",
        #"number of top links",
        #"number of top nodes",
        #"weight",
        #"weighted degree u",
        #"weighted degree v"
        "degree absolute difference",

        ]

columns0_bis = [
        #"density",
        #"average degree",
        #"average weight",
        #"average weighted degree",
        #"weighted degree absolute difference",
        ]

columns1 = ["number of links",
        "number of nodes",
        "count nodes degree 1",
        "count nodes degree 2",
        "degree absolute difference",
        "degree u",
        "degree v",
        "max degree",
        "max weighted degree",
        "number of top links",
        "number of top nodes",
        "total weight",
        "weight",
        "weighted degree absolute difference",
        "weighted degree u",
        "weighted degree v",
        "average degree",
        "average weight",
        "average weighted degree",
        "density",

        "top degree",
        "top degree v",
        "top max weighted degree",
        "top weighted degree"]
columns2 = ["Jaccard index bipartite u",
        "Jaccard index bipartite v",
        "adamic adar",
        "adamic adar bipartite u",
        "adamic adar bipartite v",

        "egonet (NNu u Nv) n (Nu u NNv) maxsize",
        "egonet (NNu u Nv) n (Nu u NNv) number of links",
        "egonet (NNu u Nv) n (Nu u NNv) number of nodes",
        "egonet (NNu u Nv) u (Nu u NNv) maxsize",
        "egonet (NNu u Nv) u (Nu u NNv) number of links",
        "egonet (NNu u Nv) u (Nu u NNv) number of nodes",
        "egonet Nu number of links",
        "egonet Nu number of nodes",
        "egonet Nu u NNu maxsize",
        "egonet Nu u NNu number of links",
        "egonet Nu u NNu number of nodes",
        "egonet Nu u Nv maxsize",
        "egonet Nu u Nv number of links",
        "egonet Nu u Nv number of nodes",
        "egonet Nv number of links",
        "egonet Nv number of nodes",
        "egonet Nv u NNv maxsize",
        "egonet Nv u NNv number of links",
        "egonet Nv u NNv number of nodes"]
columns3 = ["link component number of links",
        "link component number of nodes",
        "max distance change u",
        "max distance change v",
        "max subset size in G 0",
        "max subset size in G 1",
        "max subset size in G 2",
        "max subset size in G 3",
        "max subset size in G 4",
        "max subset size in G- 0",
        "max subset size in G- 1",
        "max subset size in G- 2",
        "max subset size in G- 3",
        "max subset size in G- 4",
        "neighborhood overlap",
        "closeness u G",
        "closeness u G-",
        "closeness v G",
        "closeness v G-",
        "distance u v in G-",
        "eccentricity u in G",
        "eccentricity u in G-",
        "eccentricity v in G",
        "eccentricity v in G-",

        "number of components G",
        "number of components G-",
        "number of distance change u",
        "number of distance change v",
        "number of nodes changing partition 0",
        "number of nodes changing partition 1",
        "number of nodes changing partition 2",
        "number of nodes changing partition 3",
        "number of nodes changing partition 4",
        "number of subsets in G 0",
        "number of subsets in G 1",
        "number of subsets in G 2",
        "number of subsets in G 3",
        "number of subsets in G 4",
        "number of subsets in G- 0",
        "number of subsets in G- 1",
        "number of subsets in G- 2",
        "number of subsets in G- 3",
        "number of subsets in G- 4",
        "size largest component G",
        "size largest component G-",
        "sum difference of distances to u and v in G",
        "sum difference of distances to u and v in G-",
        "total distance change u",
        "total distance change v",
        "u subset size in G 0",
        "u subset size in G 1",
        "u subset size in G 2",
        "u subset size in G 3",
        "u subset size in G 4",
        "u subset size in G- 0",
        "u subset size in G- 1",
        "u subset size in G- 2",
        "u subset size in G- 3",
        "u subset size in G- 4",
        "u v in same community in G 0",
        "u v in same community in G 1"
        ,"u v in same community in G 2",
        "u v in same community in G 3",
        "u v in same community in G 4",
        "u v in same community in G- 0",
        "u v in same community in G- 1",
        "u v in same community in G- 2",
        "u v in same community in G- 3",
        "u v in same community in G- 4",
        "v  subset size in G 0",
        "v  subset size in G 1",
        "v  subset size in G 2",
        "v  subset size in G 3",
        "v  subset size in G 4",
        "v  subset size in G- 0",
        "v  subset size in G- 1",
        "v  subset size in G- 2",
        "v  subset size in G- 3",
        "v  subset size in G- 4",
        "core number u in G",
        "core number u in G-",
        "core number v in G",
        "core number v in G-",
        "degeneracy G",
        "degeneracy G-",
        "number common neighbors expected from degree",
        "pagerank max G",
        "pagerank max G-"]

#columns = columns1  + columns2 + columns3

def make_split(labels, K):
    # make 10-fold cross val 
    N_fraud = 0
    N_nonFraud = 0

    all_fraud = []
    all_nonFraud = []
    num2fraud = dict()

    print('getting number of elements')
    with open(labels, 'r') as fin:
        line = next(fin)
        line = next(fin) # skip header
        while line:
            num, is_fraud = line.strip().split(',')
            #if int(num) < loLIMIT:
            #    line = next(fin)
            #    continue
            if int(num) > LIMIT:
                break
            num2fraud[int(num)] = int(is_fraud)
            if is_fraud == "1":
                N_fraud += 1
                all_fraud.append(int(num))
            else:
                N_nonFraud += 1
                all_nonFraud.append(int(num))
            try:
               line = next(fin)
            except StopIteration:
                break
    
    # downsample non frauds to have balanced classes and reduce data size
    downsampled_nonFraud = list(np.random.choice(all_nonFraud, size=N_fraud, replace=False))

    # random 10-split
    #K = 10 # TODO pass as argument to generalize ?
    splits = []
    from collections import defaultdict
    split2Nfraud = defaultdict(int)
    for split in range(K):
        splits.append([])

    for sample in downsampled_nonFraud + all_fraud:
        rdm_split = np.random.randint(low = 0, high = K)
        splits[rdm_split].append(sample)
        split2Nfraud[rdm_split] += num2fraud[sample]

    for split in range(K):
        print(f"split {split} has {split2Nfraud[split]} frauds")

    #for split in range(K):
    #    splits[split] = sorted(splits[split])

    return splits, N_fraud, N_nonFraud
    #train_index = []
    #train_fraud = []
    #train_nonFraud = []
    #test_index = []

    #print('creating train and test sets')
    #with open(labels, 'r') as fin:
    #    line = next(fin)
    #    line = next(fin)
    #    #for line in labels:
    #    while line:
    #        num, is_fraud = line.strip().split(',')
    #        if int(num) < 501:
    #            # skip first 500 samples that don't have features
    #            line = next(fin)
    #            continue

    #        #if int(num) > 100000: ... pourquoi j'ai fait ça ?!?!
    #        #    break
    #        # check if fraud or not 
    #        if is_fraud == "1":
    #            # 50%/50% split
    #            if np.random.randint(0, high=4) < 2:
    #                train_fraud.append(int(num))
    #            else:
    #                test_index.append(int(num))
    #        else:
    #            ## Downsample non frauds, then 50%/50% split
    #            if np.random.randint(0, high=4) < 2:
    #                train_nonFraud.append(int(num))
    #            else:
    #                test_index.append(int(num))
    #        try:
    #           line = next(fin)
    #        except StopIteration:
    #            break

    #downsampled_nonFraud = list(np.random.choice(train_nonFraud, size=len(train_fraud), replace=False))
    #train_index = sorted([*train_fraud , *downsampled_nonFraud])
    #return train_index, test_index, N_fraud, N_nonFraud


def train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud):

    # Both train_index and test_index should be sorted :)
    #rf =  RandomForestClassifier(max_depth=200, n_estimators=1000, n_jobs=1)
    rf =  RandomForestClassifier(n_estimators=NE, n_jobs=1, verbose=1, ccp_alpha=ccp_alpha)

    is_fraud = np.zeros(n_batch)

    print(f'learning')
    new_batch = False
    with ExitStack() as stack:
        files = [stack.enter_context(open(os.path.join(dataset, filename)))
                    for filename in csv_list]
        label = stack.enter_context(open(labels))
        lab = next(label)
        lab = next(label) # skip header
        #read_files = [fin.readlines() for fin in files]
        # skip header
        idx = 0
        bidx= 0
        train_idx = 0
        batch_dim = 0
        first_batch = True
        N_feat = 0
        #metrics = []

        for num in range(0, LIMIT): #N_fraud + N_nonFraud - 2):
            try:
                all_lab = next(label).strip().split(',')
            except StopIteration:
                continue
            cur_idx = 0

            lab = all_lab[1]
            #if num <=498:
            if num <= window_max -1:
                #if num < loLIMIT:
                # skip first 500 samples, they don't have features
                continue
                 
            for num_file, fin in enumerate(files):
                if num == window_max:
                    metrics = []
                    line = next(fin).strip().split(',') # skip header
                    interaction_id = line.index('interaction_id')
                    metrics.append(interaction_id)
                    batch_dim += 1
                    #batch_dim += len(line) - 1
                    # get m1 metrics
                    for m1 in columns:
                        try:
                            metrics.append(line.index(m1)) #get metric index
                            batch_dim += 1
                        except ValueError:
                            #print(f'did not find {m1}')
                            pass

                try:
                    _features = next(fin).strip().split(',')
                    __features = [_features[m1] for m1 in metrics]
                    #_features = _features[:14] + _features[33:]

                except StopIteration:
                    break
                #it_id = int(_features[interaction_id])
                it_id = int(__features[0])
                
                # take greatest window size as first interaction ID for ALL windows
                while (it_id < window_max):
                    try:
                        _features = next(fin).strip().split(',')
                        __features = [_features[m1] for m1 in metrics]
                        #_features = _features[:14] + _features[33:]

                    except StopIteration:
                        break
                    #it_id = int(_features[interaction_id])
                    it_id = int(__features[0])


                #if it_id != int(all_lab[0]) and it_id != int(all_lab[0]) + 1:
                #    print(f'{it_id} {all_lab[0]}')
                #if  it_id != int(all_lab[0]):
                while (it_id > int(all_lab[0])):
                    # Handle cases where we skipped a n interaction because it was a loop
                    try:
                        all_lab = next(label).strip().split(',')
                    except StopIteration:
                        break
                if it_id != int(all_lab[0]):
                    print(f'{it_id} {all_lab[0]}')

                features = np.array([float(f) for f in __features])
                while (train_idx < len(train_index) and it_id > train_index[train_idx]):
                    train_idx += 1
                if not first_batch and train_idx < len(train_index) and it_id == train_index[train_idx]:
                    #else:
                    N_feat = len(features)
                    batch[idx % n_batch, cur_idx : cur_idx + N_feat ] = features
                    cur_idx += N_feat

                #else:
                #    continue
                if num == window_max and first_batch:
                    batch = np.zeros((n_batch, batch_dim))
                    first_batch = False
                    continue
                if train_idx < len(train_index) and  it_id == train_index[train_idx]:
                    is_fraud[idx % n_batch] = lab
                    idx += 1
                    #train_index.pop(0)
                    train_idx += 1
                    new_batch = True
                # when batch is ready, learn !
                if idx  % n_batch  == 0 and new_batch == True:
                    bidx += 1
                    rf.fit(batch, is_fraud)
                    #fraud_pred = rf.predict(batch)
                    #roc = roc_auc_score(fraud_pred, is_fraud)

                    print(f'training batch score {rf.score(batch, is_fraud)}')

                    new_batch = False # needed because next line might not be in train_index
        else:
            if idx % n_batch > 0:
                rf.fit(batch[:idx % n_batch,:], is_fraud[:idx % n_batch] )
                fraud_pred = rf.predict(batch[:idx % n_batch,:])

                print(f'training batch score {rf.score(batch, is_fraud)}')
    return rf

def predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf):

    # predict
    print("predicting")
    #batch = np.zeros((n_batch, N_feat * len(csv_list)))
    bidx = 0
    y_test = []
    new_batch = False
    with ExitStack() as stack:
        files = [stack.enter_context(open(os.path.join(dataset, filename)))
                    for filename in csv_list]
        idx = 0
        label = stack.enter_context(open(labels))
        lab = next(label)
        lab = next(label) # skip header
        test_idx = 0
        first_batch = True
        batch_dim = 0
        #metrics = []

        # skip header
        for num in range(0, LIMIT):# N_fraud + N_nonFraud-1):
            try:
                all_lab = next(label).strip().split(',')
            except StopIteration:
                continue
            cur_idx = 0

            lab = all_lab[1]

            if num <window_max :
                continue
            for num_file, fin in enumerate(files):
                if num == window_max:
                    metrics = []
                    line = next(fin).strip().split(',') # skip header
                    interaction_id = line.index('interaction_id')
                    #batch_dim += len(line) - 1
                    metrics.append(interaction_id)
                    #batch_dim += len(line) - 1
                    batch_dim += 1
                    # get m1 metrics
                    for m1 in columns:
                        try:
                            metrics.append(line.index(m1)) #get metric index
                            batch_dim += 1
                        except ValueError:
                            #print(f'did not find {m1}')
                            pass

                try:
                    _features = next(fin).strip().split(',')
                    __features = [_features[m1] for m1 in metrics]

                    #_features = _features[:14] + _features[33:]
                except StopIteration:
                    break
                #it_id = int(_features[interaction_id] )
                it_id = int(__features[0])

                #if  it_id != int(all_lab[0]):
                while (it_id > int(all_lab[0])):
                    # Handle cases where we skipped a n interaction because it was a loop
                    try:
                        all_lab = next(label).strip().split(',')
                    except StopIteration:
                        break
                if it_id != int(all_lab[0]):
                    print(f'{it_id} {all_lab[0]}')
                features = np.array([float(f) for f in __features])
                while (test_idx < len(test_index) and it_id > test_index[test_idx]):
                    test_idx += 1
                if not first_batch and test_idx < len(test_index) and  it_id == test_index[test_idx]:
                    #if first_batch:
                    #    batch = np.zeros((n_batch, batch_dim))
                    #    first_batch = False
                    #else:
                    N_feat = len(features)

                    batch[idx % n_batch, cur_idx : cur_idx + N_feat ] = np.array(features)
                    #batch[idx  % n_batch, (N_feat * num_file) : (N_feat * (num_file +1)) ] = features
                    cur_idx += N_feat

                #else:
                #    continue
                if num == window_max and first_batch:
                    batch = np.zeros((n_batch, batch_dim))
                    first_batch = False
                    continue
                if test_idx < len(test_index) and  it_id == test_index[test_idx]:
                    y_test.append(int(lab))
                    idx += 1
                    #test_index.pop(0)
                    test_idx += 1
                    new_batch = True
                # when batch is ready, learn !
                if idx  % n_batch  == 0 and new_batch == True:
                    _y_pred = rf.predict(batch)
                    if bidx == 0:
                        y_pred = _y_pred
                    else:
                        y_pred = np.concatenate([y_pred, _y_pred])
                    bidx += 1
                    new_batch = False
        else:
            if idx % n_batch > 0:

                _y_pred = rf.predict(batch[:idx % n_batch,:])
                if bidx == 0:
                    y_pred = _y_pred
                else:
                    y_pred = np.concatenate([y_pred, _y_pred])
    return y_pred, y_test

def score(y_pred, y_test):
    roc = roc_auc_score(np.array(y_test), y_pred)
    recall =recall_score(y_test, y_pred)
    accuracy =accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    f1 = f1_score(y_test,y_pred)
    print('forest  recall {} accuracy {} precision {} f1 {} roc {}'.format( recall, accuracy, precision, f1, roc))
    return roc, recall, accuracy, precision, f1

def by_columns(csv_list, dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    global columns
    col2grp = {0: "O(1)", 1:"O(k)", 2:"O(n+m)"}

    #train_index, test_index, N_fraud, N_nonFraud = split(args.labels)

    #print(f'{train_index[0]} {test_index[0]} {len(train_index)} {len(test_index)}')
    #csv_list = []
    #csv_list = _csv_list
    #csv_list_unit = _csv_list
    #csv = csv_list[-1]
    col_cum = []
    for col_idx, col in enumerate([columns1, columns2, columns3]):

        #csv_list.append(csv)
        #csv_list_unit = [csv]
        col_cum += col 
        #col_cum.append(col)

        split_scores_unit = []
        split_scores_all = []

        #print(col2grp[col_idx])
        for split_idx, split in enumerate(splits):
            test_index = sorted(split)
            train_index = []
            for k in range(len(splits)):
                if k == split_idx:
                    continue
                train_index += splits[k]
            train_index = sorted(train_index)
            columns = col
            rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)

            y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {col2grp[col_idx]}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_unit.append((roc,recall, accuracy, precision, f1))


            columns = col_cum
            rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)
            y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {" ".join(col_cum)}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_all.append((roc,recall, accuracy, precision, f1))
        #with open(os.path.join(output , f'{col2grp[col_idx]}_kinf_allwindows.txt'), 'w' )as fout:
        #    fout.write(f'roc,recall,accuracy,precision,f1\n')
        #    for roc, recall, acc, prec, f1 in split_scores_unit:
        #        fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

        #with open(os.path.join(output + f'all_{col2grp[col_idx]}_kinf_allwindows.txt'), 'w') as fout:
        #    fout.write(f'roc,recall,accuracy,precision,f1\n')
        #    for roc, recall, acc, prec, f1 in split_scores_all:
        #        fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

def by_windows(_csv_list, dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    global columns
    columns = columns1 + columns2 + columns3
    #col2grp = {0: "O(1)", 1:"O(k + k*k)", 2:"O(n+m) and others")
    #train_index, test_index, N_fraud, N_nonFraud = split(args.labels)

    #print(f'{train_index[0]} {test_index[0]} {len(train_index)} {len(test_index)}')
    csv_list = []
    csv_list_unit = []
    #csv_list = _csv_list
    #csv_list_unit = _csv_list
    #csv = csv_list[-1]
    #for col_idx, col in enumerate([columns1, columns2, columns3]):
    for csv in _csv_list:
        
        csv_list.append(csv)
        csv_list_unit = [csv]
        split_scores_unit = []
    
        split_scores_cum = []

        for split_idx, split in enumerate(splits):
            test_index = sorted(split)
            train_index = []
            for k in range(len(splits)):
                if k == split_idx:
                    continue
                train_index += splits[k]
            train_index = sorted(train_index)
            rf = train(dataset, labels, csv_list_unit, n_batch, train_index,  N_fraud, N_nonFraud)

            y_pred, y_test = predict(dataset, labels, csv_list_unit, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {csv}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_unit.append((roc,recall, accuracy, precision, f1))

            rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)
            y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {".".join(csv_list)}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_cum.append((roc,recall, accuracy, precision, f1))

        with open(os.path.join(output, f'{csv}_kinf_allcolumns.txt'), 'w') as fout:
            fout.write(f'roc,recall,accuracy,precision,f1\n')
            for roc, recall, acc, prec, f1 in split_scores_unit:
                fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

        with open(os.path.join(output, f'{".".join(csv_list)}_kinf_allcolumns.txt'), 'w') as fout:
            fout.write(f'roc,recall,accuracy,precision,f1\n')
            for roc, recall, acc, prec, f1 in split_scores_cum:
                fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

def by_deg_bound(csv_list, _dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    global columns
    print('by_deg')
    #columns = columns1 + columns2 + columns3
    columns = columns2
    #for col_idx, col in enumerate([columns1, columns2, columns3]):
    
    #for csv in _csv_list:
    for k_bound in [1, 10, 100, 1000, 100000]:
        dataset = os.path.join(_dataset, f"k{k_bound}")
        split_scores = []

        for split_idx, split in enumerate(splits):
            test_index = sorted(split)
            train_index = []
            for k in range(len(splits)):
                if k == split_idx:
                    continue
                train_index += splits[k]
            train_index = sorted(train_index)
            rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)

            y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {k_bound}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores.append((roc, recall, accuracy, precision, f1))
        with open(os.path.join(output, f'k{k_bound}_allwindows_allcolumns.txt'), 'w') as fout:
            fout.write(f'roc,recall,accuracy,precision,f1\n')
            for roc, recall, acc, prec, f1 in split_scores:
                fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

def complete(csv_list, dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    global columns
    col2grp = {0: "O(1)", 1:"O(k)", 2:"O(n+m)"}

    #train_index, test_index, N_fraud, N_nonFraud = split(args.labels)

    #print(f'{train_index[0]} {test_index[0]} {len(train_index)} {len(test_index)}')
    #csv_list = []
    #csv_list = _csv_list
    #csv_list_unit = _csv_list
    #csv = csv_list[-1]
    col_cum = []
    #for col_idx, col in enumerate([columns1, columns2, columns3]):
    columns = columns1 + columns2 + columns3

    #csv_list.append(csv)
    #csv_list_unit = [csv]
    #col_cum += col 
    #col_cum.append(col)

    split_scores_unit = []
    split_scores_all = []

    #print(col2grp[col_idx])
    names = ['interaction_id'] + columns
    for split_idx, split in enumerate(splits):
        test_index = sorted(split)
        train_index = []
        for k in range(len(splits)):
            if k == split_idx:
                continue
            train_index += splits[k]
        train_index = sorted(train_index)
        rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)

        # print interpretation
        #interpretation = sorted(zip(map(lambda x: round(x, 6), rf.feature_importances_), names), reverse=True)
        #sum_scores = np.sum(rf.feature_importances_)
        #print(f'- sum importance {sum_scores}')
        #for key, val in interpretation:
        #    print(f'- {val} {key}')

        y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
        roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
        split_scores_unit.append((roc,recall, accuracy, precision, f1))

    #with open(os.path.join(output , f'{col2grp[col_idx]}_kinf_allwindows.txt'), 'w' )as fout:
    #    fout.write(f'roc,recall,accuracy,precision,f1\n')
    #    for roc, recall, acc, prec, f1 in split_scores_unit:
    #        fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

    #with open(os.path.join(output + f'all_{col2grp[col_idx]}_kinf_allwindows.txt'), 'w') as fout:
    #    fout.write(f'roc,recall,accuracy,precision,f1\n')
    #    for roc, recall, acc, prec, f1 in split_scores_all:
    #        fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", help="path to dataset")
    parser.add_argument("labels", help="path to labels")
    parser.add_argument("output", help="path to output")
    parser.add_argument("-k",'--K', type=int, default=10, help="number of splits")
    group = parser.add_mutually_exclusive_group()

    group.add_argument("-c", "--by_columns", default=False, action="store_true", help="compare metrics by columns")
    group.add_argument("-w", "--by_windows", default=False, action="store_true", help="compare metrics by window sizes")
    group.add_argument("-bk", "--by_limit",  default=False, action="store_true", help="compare metrics by limit k")
    group.add_argument('--complete', default=False, action="store_true", help="all parameters to the max")
    
    #parser.add_argument("N_feat", type=int, help="number of features")
    #parser.add_argument("interaction_id", type=int, help="number of features")
    #from IPython.core.debugger import Tracer
    #Tracer()() #this one triggers the debugger

    args = parser.parse_args()

    mySeed = np.random.randint(2**32-1)
    #mySeed= 2659021839
    np.random.seed(mySeed) # set random seed
    print(f'numpy random seed set at {mySeed}')
    n_batch = 40000

    # list of csvs to read for features
    csv_list = [
    #"Output_H_8.csv",
    #"Output_H_64.csv",
    #"Output_H_128.csv",
    #"Output_H_512.csv",
    #"Output_H_1024.csv",
    #"Output_H_2048.csv",
    #"Output_H_4096.csv",
    "Output_H_8192.csv",#,
    #"Output_H_65536.csv",
    #"Output_H_65636.csv",
    #"Output_H_131072.csv",#]#,
    #"Output_H_1048576.csv"#]
    ]
    global window_max
    windows = [int( csv.split('_')[2].split('.')[0] ) for csv in csv_list]
    window_max = max(windows)
    #csv_list = ['Output_H_512.csv','Output_H_4096.csv', 'Output_H_131072.csv']
    global NE 
    NE = 1000
    global ccp_alpha
    global LIMIT
    global loLIMIT
    #LIMIT = 200000
    LIMIT = 100000000
    ccp_alpha=0.03
    splits, N_fraud, N_nonFraud = make_split(args.labels, args.K)
    if args.by_columns:
        by_columns(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)
    elif args.by_windows:
        by_windows(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)
    elif args.by_limit:
        by_deg_bound(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)
    elif args.complete:
        splits, N_fraud, N_nonFraud = make_split(args.labels, args.K)
        #LIMIT = N_fraud + N_nonFraud
        #for NE in [100, 1000, 10000, 100000]:
        print(NE)
        #for LIMIT in [10000, 50000, 100000, 250000, 500000, 750000, 1000000]:
        #for ccp_alpha in np.linspace(0,0.1,11):
        for ccp_alpha in [0.03]:
            print(f'ccp_alpha {ccp_alpha}')
            complete(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)

if __name__ == "__main__":
    main()
