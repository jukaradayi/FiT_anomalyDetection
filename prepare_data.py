import os
import sys
import wget
import gzip
import shutil
import random
import datetime
import numpy as np
import pandas as pd
from subprocess import check_call


"""

## Usage

To download the primaryschool dataset from Socio-Pattern and format it
into a proper stream with integer node labels, run:

```
$ python3 prepare_data.py primaryschool primaryschool mapping.txt
```

This downloads the data (./primaryschool.csv.gz), 
extract the data (./primaryschool.csv),
and build a proper tuv stream (./primaryschool.gz),
with node mapping (./mapping.txt) 

"""


def get_primary_school_data(directory):
    """
    Retrieve the Primary school dataset from sociopatterns.
    """
    url = "http://www.sociopatterns.org/wp-content/uploads/2015/09/primaryschool.csv.gz"
    print("Download the data set from {url}".format(url=url))
    filename = wget.download(url, directory)
    with gzip.open(filename, 'rb') as f_in:
        with open(os.path.join(directory, "primaryschool.csv"), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    # Build the dataframe...
    df = pd.read_csv(os.path.join(directory, "primaryschool.csv"), header=None, sep="\t")
    return df

def build_columns(df, use_datetime=False):
    """
    Build the dataframe columns.
    """
    df['Source'] = df[1].astype(str)
    df['Target'] = df[2].astype(str)
    if use_datetime:
        # Set a fake starting date since the data starts at zero...
        original_date = pd.to_datetime("2020-01-01 00:00:00")
        daily_date = original_date
        df["Datetime"] = pd.to_datetime(df[0].apply(
                                lambda x: daily_date + datetime.timedelta(seconds=x)).apply(
                                lambda x: x.strftime('%Y-%m-%d %H:%M:%S')))
    else:
        df['Datetime'] = df[0].astype(int)
    df = df[['Datetime', 'Source', 'Target']]
    return df


def load_data(dataset):
    """Load or download data."""
    if os.path.basename(dataset) == "primary_school.csv":
        #stream_filename = "./primary_school.csv"

        # If we already have the dataset in the current folder, then use it
        if os.path.exists(dataset):
            print(f"Using previously downloaded dataset: {stream_filename}")
            df = pd.read_csv(dataset, header=None, sep="\t")
        # Otherwise download it.    
        else:
            print(f"Downloading dataset and writing to disk as {stream_filename}")
            df = get_primary_school_data(".")
        df = build_columns(df)
        return df
    else:
        # Assume (t,u,v) with t in seconds, no header, separated by tabs
        return pd.read_csv(dataset, header=None, sep=" ", names=['Datetime', 'Source', 'Target'])


def write_stream(stream, path, mapping_path):
    """
    Write the given stream at the given location.
    """
    if os.path.exists(path):
        print(f"WARNING: File {path} already exists. Overwrite...")
    print(f"Writing stream to {path}...")
    last_idx = 0
    mapping = dict()
    with open(path, "w") as fp:
        for idx,interaction in stream.iterrows():
            t,u,v = interaction['Datetime'], interaction['Source'], interaction['Target']
            idxs = [-1, -1]            
            if idx==0:
                t0 = t
            for i,x in enumerate([u,v]):
                if x in mapping:
                    idxs[i] = mapping[x]
                else:
                    mapping[x] = last_idx
                    idxs[i] = last_idx
                    last_idx += 1
            fp.write(f"{int((t-t0))} {idxs[0]} {idxs[1]}\n")
    with open(mapping_path, "w") as fp:
        for k,v in mapping.items():
            fp.write(f"{k} {v}\n")

def main(dataset, output_stream, output_mapping):
    df = load_data(dataset)
    write_stream(df, output_stream, output_mapping)
    #check_call(['gzip', output_stream])
    return 0

if __name__ == "__main__":
    dataset = sys.argv[1]
    output_stream = sys.argv[2]
    output_mapping = sys.argv[3]
    main(dataset, output_stream, output_mapping)
