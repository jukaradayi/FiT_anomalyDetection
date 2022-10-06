import os
import sys
import ipdb
import argparse
import numpy as np

from collections import defaultdict
from contextlib import ExitStack
from sklearn.ensemble import RandomForestClassifier #, AdaBoostClassifier
from sklearn.metrics import f1_score, recall_score, accuracy_score, precision_score, roc_auc_score, mean_squared_error, roc_curve, plot_roc_curve

# Categories of metrics -- Useful to compare influence of group of columns on results
columns0 = [
        "number of links",
        "number of nodes",
        "total weight",
        "count nodes degree 1",
        "count nodes degree 2",
        "degree u",
        "degree v",

        "max degree",
        "max weighted degree",
        "number of top links",
        "number of top nodes",
        "weight",
        "weighted degree u",
        "weighted degree v"
        "degree absolute difference",

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


def make_split(labels, K):
    ''' Create k splits for k-fold cross validation.
        Loop once through labels to get samples of "Fraud" and "Non-Fraud",
        then downsample the list of "Non-Fraud" to have the same number as "Fraud",
        then place each sample of "Fraud" and "Non-Fraud" randomly into one of the k splits.
    '''
    # Number of Fraud and Non-Fraud
    N_fraud = 0
    N_nonFraud = 0

    # List of fraud and non-fraud
    all_fraud = []
    all_nonFraud = []
    num2fraud = dict() # mapping from interaction-ID to label

    print('getting number of elements')
    with open(labels, 'r') as fin:
        line = next(fin)
        line = next(fin) # skip header
        while line:
            num, is_fraud = line.strip().split(',')
            ## loLimit = 0 by default, can be used to train only on subset of dataset 
            #if int(num) < loLIMIT:
            #    line = next(fin)

            #    continue

            # window_max is the size of the largest window used to extract features
            # skip interaction_ids < window_max to use only interaction_id where 
            # all window have the feature extracted.
            if int(num) < window_max:
                line = next(fin)
                continue

            # By default LIMIT > len(fin), can be used to train only on subset of dataset
            #if int(num) > LIMIT: 
            #    break

            # map interaction_id to label, and count Frauds and Non-Frauds
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

    splits = []
    split2Nfraud = defaultdict(int) # map split ID to its number of Frauds

    for split in range(K):
        splits.append([])

    for sample in downsampled_nonFraud + all_fraud:
        # randomly pick split for each sample
        rdm_split = np.random.randint(low = 0, high = K)
        splits[rdm_split].append(sample)
        split2Nfraud[rdm_split] += num2fraud[sample]

    # print number of fraud by split
    for split in range(K):
        print(f"split {split} has {split2Nfraud[split]} frauds")

    return splits, N_fraud, N_nonFraud


def train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud):
    ''' Train random forest on the data.
    
    Arguments:
    ----------
    dataset: str, path to the csv files containing the features
    labels: str, path to the csv file containing the labels 
    csv_list: list(str), name of the csv used for training
    n_batch: int, size of the batches to feed to the random forest 
    train_index: list(int), list of the interaction_ids to use for this training
    N_fraud: int, number of fraud
    N_nonFraud: int, number of non-fraud

    Return:
    -------
    rf: trained random forest
    '''
    # Both train_index and test_index should be sorted :)
    # number of estimators, max depth of the tree and pruning coefficient ccp_alpha are defined in main()
    # as global variables
    rf =  RandomForestClassifier(n_estimators=NE, n_jobs=1, verbose=0, max_depth=MD,  ccp_alpha=ccp_alpha)

    is_fraud = np.zeros(n_batch)

    print(f'learning')
    new_batch = False

    # start training
    with ExitStack() as stack:

        # read all csvs in parallel
        files = [stack.enter_context(open(os.path.join(dataset, filename)))
                    for filename in csv_list]
        label = stack.enter_context(open(labels))
        lab = next(label)
        lab = next(label) # skip header


        idx = 0 # line number
        bidx= 0 # batch number
        train_idx = 0 # index of train_index list
        batch_dim = 0 # number of features used (sum of N_feat for each csv file)
        first_batch = True # boolean to indicate if it is the first_batch
        N_feat = 0 # number of features for current csv 
        metrics = defaultdict(list) # for each csv file keep the list of metrics

        for num in range(0, LIMIT):
            try:
                all_lab = next(label).strip().split(',')
            except StopIteration:
                continue
            cur_idx = 0

            lab = all_lab[1]

            # skip lines for which we don't have all features
            if num <= window_max :
                continue

            # loop through each csv file to fill current batch.
            for num_file, fin in enumerate(files):
                if num == window_max +1  :

                    line = next(fin).strip().split(',') # skip header
                    interaction_id = line.index('interaction_id')
                    metric_it_id = interaction_id

                    # get the position of each metric and count increment dimension of batch
                    for m1 in columns:
                        try:
                            metrics[num_file].append(line.index(m1)) #get metric index
                            batch_dim += 1
                        except ValueError:
                            # might happen if a feature is not present in one of the CSV
                            pass

                # read feature in current csv
                try:
                    _features = next(fin).strip().split(',')
                    __features = [_features[m1] for m1 in metrics[num_file]]
                except StopIteration: # StopIteration exception when reaching end of csv file
                    break

                it_id = int(_features[metric_it_id]) # get current interaction_id
                
                # take greatest window size as first interaction ID for ALL windows
                while (it_id < int(all_lab[0])):
                    try:
                        _features = next(fin).strip().split(',')
                        __features = [_features[m1] for m1 in metrics[num_file]]
                    except StopIteration:
                        break
                    it_id = int(_features[metric_it_id])

                # Handle cases where we skipped a n interaction because it was a loop
                while (it_id > int(all_lab[0])):
                    try:
                        all_lab = next(label).strip().split(',')
                    except StopIteration:
                        break

                if it_id != int(all_lab[0]):
                    # interaction_id don't match between labels and features
                    # Shouldn't happen
                    print(f'WARNING: mismatch, feature: {it_id} label: {all_lab[0]}')

                features = np.array([float(f) for f in __features])
                while (train_idx < len(train_index) and it_id > train_index[train_idx]):
                    train_idx += 1
                
                # fill batch
                if not first_batch and train_idx < len(train_index) and it_id == train_index[train_idx]:
                    N_feat = len(features)
                    batch[idx % n_batch, cur_idx : cur_idx + N_feat ] = features
                    cur_idx += N_feat

                # initialize first batch
                if num == window_max +1 and first_batch:
                    batch = np.zeros((n_batch, batch_dim))
                    first_batch = False
                    continue

                if train_idx < len(train_index) and  it_id == train_index[train_idx]:
                    is_fraud[idx % n_batch] = lab
                    idx += 1
                    train_idx += 1
                    new_batch = True

                # when batch is ready, learn !
                if idx  % n_batch  == 0 and new_batch == True:
                    bidx += 1
                    rf.fit(batch, is_fraud)
                    # print training score for current batch
                    print(f'training batch score {rf.score(batch, is_fraud)}')
                    new_batch = False # needed because next line might not be in train_index
        else:
            # train on last batch
            if idx % n_batch > 0:
                rf.fit(batch[:idx % n_batch,:], is_fraud[:idx % n_batch] )
                fraud_pred = rf.predict(batch[:idx % n_batch,:])
                print(f'training batch score {rf.score(batch, is_fraud)}')

    return rf

def predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf):
    ''' Launch prediction with trained random forest.
        See train for argument list.
    '''
    # predict
    print("predicting")
    bidx = 0 # batch index
    y_test = [] # labels
    interactions = [] # write interaction ids into file - Useful to study features of frauds that were correctly detected
    new_batch = False

    # Same as in training...
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
        metrics = defaultdict(list)

        # skip header
        for num in range(0, LIMIT):

            try:
                all_lab = next(label).strip().split(',')
            except StopIteration:
                continue
            cur_idx = 0
            lab = all_lab[1]

            if num <= window_max :
                continue

            for num_file, fin in enumerate(files):
                if num == window_max+1:
                    line = next(fin).strip().split(',') # skip header
                    interaction_id = line.index('interaction_id')
                    metric_it_id = interaction_id

                    # get metrics
                    for m1 in columns:
                        try:
                            metrics[num_file].append(line.index(m1)) #get metric index
                            batch_dim += 1
                        except ValueError:
                            pass

                try:
                    _features = next(fin).strip().split(',')
                    __features = [_features[m1] for m1 in metrics[num_file]]
                except StopIteration:
                    break
                it_id = int(_features[metric_it_id])

                # take greatest window size as first interaction ID for ALL windows
                while (it_id < int(all_lab[0])):
                    try:
                        _features = next(fin).strip().split(',')
                        __features = [_features[m1] for m1 in metrics[num_file]]

                    except StopIteration:
                        break
                    it_id = int(_features[metric_it_id])



                while (it_id > int(all_lab[0])):
                    # Handle cases where we skipped a n interaction because it was a loop
                    try:
                        all_lab = next(label).strip().split(',')
                    except StopIteration:
                        break

                if it_id != int(all_lab[0]):
                    # shouldn't happen
                    print(f'WARNING: mismatch, feature: {it_id} label: {all_lab[0]}')

                features = np.array([float(f) for f in __features])
                while (test_idx < len(test_index) and it_id > test_index[test_idx]):
                    test_idx += 1

                if not first_batch and test_idx < len(test_index) and  it_id == test_index[test_idx]:
                    N_feat = len(features)
                    batch[idx % n_batch, cur_idx : cur_idx + N_feat ] = np.array(features)
                    cur_idx += N_feat

                if num == window_max +1 and first_batch:
                    batch = np.zeros((n_batch, batch_dim))
                    first_batch = False
                    continue

                if test_idx < len(test_index) and  it_id == test_index[test_idx]:
                    y_test.append(int(lab))
                    interactions.append(it_id)
                    idx += 1
                    test_idx += 1
                    new_batch = True

                # when batch is ready, PREDICT !
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

    # Useful for studying features of correctly predicted fraud
    #with open('bleckwen_predictions', 'w') as fout:

    #    for pred, truth, idx in zip(y_pred, y_test, interactions):
    #        if truth == 1 and pred == 1:
    #            fout.write(f'{idx}\n')

    return y_pred, y_test

def score(y_pred, y_test):
    ''' Compute ROC-AUC score on prediction.
    Arguments:
    ---------
    y_pred: list(int): predictions, 1 if fraud, else 0
    y_test: list(int): labels, 1 if fraud, else 0
    '''
    # ROC-AUC
    roc = roc_auc_score(np.array(y_test), y_pred)

    # recall, accuracy, precision, f1 score...
    recall =recall_score(y_test, y_pred)
    accuracy =accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    f1 = f1_score(y_test,y_pred)

    print('forest  recall {} accuracy {} precision {} f1 {} roc {}'.format( recall, accuracy, precision, f1, roc))

    return roc, recall, accuracy, precision, f1

def by_columns(csv_list, dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    ''' Study the effect of each category of column in the csv file.
        Iteratively add group of column and train on these features.
    '''
    global columns
    col2grp = {0: "O(1)", 1:"O(k)", 2:"O(n+m)"}

    col_cum = [] # train cumulatively on each column
    for col_idx, col in enumerate([columns1, columns2, columns3]):
        #col_idx = 1
        #for col in [columns2]:

        col_cum += col 

        split_scores_unit = []
        split_scores_all = []

        for split_idx, split in enumerate(splits):
            test_index = sorted(split)
            train_index = []
            for k in range(len(splits)):
                if k == split_idx:
                    continue
                train_index += splits[k]
            train_index = sorted(train_index)

            # Train only on current column set
            columns = col
            rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)
            y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {col2grp[col_idx]}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_unit.append((roc,recall, accuracy, precision, f1))

            # Train on current column set AND previous column sets
            columns = col_cum
            rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)
            y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {" ".join(col_cum)}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_all.append((roc,recall, accuracy, precision, f1))

def by_windows(_csv_list, dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    ''' Study the effect of each window size on the results.
    '''
    
    global columns
    columns = columns1 + columns2 + columns3
    csv_list = []
    csv_list_unit = []
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

            # train only on current window size
            rf = train(dataset, labels, csv_list_unit, n_batch, train_index,  N_fraud, N_nonFraud)

            y_pred, y_test = predict(dataset, labels, csv_list_unit, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {csv}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_unit.append((roc,recall, accuracy, precision, f1))

            # train on current window size AND previous window sizes
            rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)
            y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
            print(f'result for {".".join(csv_list)}')
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores_cum.append((roc,recall, accuracy, precision, f1))

        #with open(os.path.join(output, f'{csv}_kinf_allcolumns.txt'), 'w') as fout:
        #    fout.write(f'roc,recall,accuracy,precision,f1\n')
        #    for roc, recall, acc, prec, f1 in split_scores_unit:
        #        fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

        #with open(os.path.join(output, f'{".".join(csv_list)}_kinf_allcolumns.txt'), 'w') as fout:
        #    fout.write(f'roc,recall,accuracy,precision,f1\n')
        #    for roc, recall, acc, prec, f1 in split_scores_cum:
        #        fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

def by_deg_bound(csv_list, _dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    ''' Study the effect of the degree bound on the results
    '''
    global columns
    columns = columns1 + columns2 + columns3
    #columns = columns2
    #for col_idx, col in enumerate([columns1, columns2, columns3]):
    
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
            roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
            split_scores.append((roc, recall, accuracy, precision, f1))
        #with open(os.path.join(output, f'k{k_bound}_allwindows_allcolumns.txt'), 'w') as fout:
        #    fout.write(f'roc,recall,accuracy,precision,f1\n')
        #    for roc, recall, acc, prec, f1 in split_scores:
        #        fout.write(f'{roc},{recall},{acc},{prec},{f1}\n')

def complete(csv_list, dataset, splits, N_fraud, N_nonFraud, labels, output, n_batch):
    '''Complete training
    '''
    global columns
    col2grp = {0: "O(1)", 1:"O(k)", 2:"O(n+m)"}

    col_cum = []
    columns = columns1 + columns2 + columns3
    split_scores_unit = []
    split_scores_all = []

    names = ['interaction_id'] + columns
    for split_idx, split in [(0,splits[0])]:
        test_index = sorted(split)
        train_index = []
        for k in range(len(splits)):
            if k == split_idx:
                continue
            train_index += splits[k]
        train_index = sorted(train_index)
        rf = train(dataset, labels, csv_list, n_batch, train_index,  N_fraud, N_nonFraud)
        y_pred, y_test = predict(dataset, labels, csv_list, n_batch, test_index, N_fraud, N_nonFraud, rf)
        roc, recall, accuracy, precision, f1 = score(y_pred, y_test)
        split_scores_unit.append((roc,recall, accuracy, precision, f1))

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
    
    args = parser.parse_args()

    mySeed = np.random.randint(2**32-1)
    np.random.seed(mySeed) # set random seed
    print(f'numpy random seed set at {mySeed}')
    n_batch = 400000

    # list of csvs to read for features
    csv_list = [
    #"Output_H_8.csv",
    #"Output_H_64.csv",
    #"Output_H_128.csv",
    #"Output_H_512.csv",
    #"Output_H_1024.csv",

    #"Output_G_1024.csv",

    ##"Output_H_2048.csv",
    ##"Output_H_4096.csv",
    "Output_H_8192.csv",#,
    #"Output_G_86400.csv",
    #"Output_H_65536.csv",
    #"Output_m3_512.csv",
    ##"Output_H_65636.csv",
    #"Output_H_131072.csv",#]#,
    #"Output_H_524288.csv",
    #"Output_H_1048576.csv"#]
    ]

    # Set global variables, in particular some for the random forest.
    global window_max
    windows = [int( csv.split('_')[2].split('.')[0] ) for csv in csv_list]
    window_max = max(windows)
    global NE # number of estimators for the random forest
    global MD # max depth of the random forest  decision trees
    #NE = 5000
    #NE = 5000
    NE = 100
    #MD = 50

    global ccp_alpha # pruning coefficient for the random forest
    global LIMIT
    global loLIMIT
    #loLIMIT = 12000000
    loLIMIT = 0
    #LIMIT = 14000000
    LIMIT = 14088142
    #LIMIT = 100000
    #ccp_alpha=0.03
    #splits, N_fraud, N_nonFraud = make_split(args.labels, args.K)
    if args.by_columns:
        ccp_alpha = 0.0001
        MD = 100
        NE = 100
        by_columns(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)
    elif args.by_windows:
        by_windows(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)
    elif args.by_limit:
        by_deg_bound(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)
    elif args.complete:
        #for MULTIPLE in [1, 2, 3, 4, 5, 7, 9]:
        splits, N_fraud, N_nonFraud = make_split(args.labels, args.K)
        #LIMIT = N_fraud + N_nonFraud
        #for NE in [100, 1000, 10000, 100000]:
        print(NE)
        #for LIMIT in [10000, 50000, 100000, 250000, 500000, 750000, 1000000]:
        #for ccp_alpha in np.linspace(0,0.1,11):
        #for ccp_alpha in [0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.03] : #, 0.03, 0.05, 0.1, 0.5]:
        for ccp_alpha in [0.0001]:
            #for MD in [25, 50, 75100, 200, 500, 1000]:
            #for MD in [None, 50, 100]: 
            for MD in [100]:
                #for NE in [100, 1000, 2000, 3000, 5000, 10000]:
                for NE in [100]:
                    #for ccp_alpha in [0.03]:   
                    print(f'MD {MD} ccp_alpha {ccp_alpha} NE {NE}')
                    complete(csv_list, args.dataset, splits, N_fraud, N_nonFraud, args.labels, args.output, n_batch)

if __name__ == "__main__":
    main()
