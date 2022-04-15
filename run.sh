#!/bin/bash
#
# run fit script with required parameters
export LD_LIBRARY_PATH=/home/fit/Workspace/networkit_lib/lib; # to get libnetworkit.so
# input arguments
set -e
export LC_ALL=C

if [ "$#" -lt 2 ]; then
    echo "Usage:"
    echo "$0 <input_file> <njobs>"
    echo "Please modify run.sh to specify window sizes"
    echo "Outputs will be written in input folder"
fi

INPUT=$1
NJOBS=$2

if [ $(($NJOBS%2)) -eq 0 ]
then
    HJOBS=$(($NJOBS / 2))
    GJOBS=$(($NJOBS / 2))
else
    HJOBS=$(($NJOBS / 2))
    GJOBS=$(($NJOBS / 2 + 1))
fi
   
echo "$HJOBS $GJOBS"

# Set H and G params
HPARAMS="8 64 512 4096 16384"
GPARAMS="10 60 3600 86400"

## set paths
FITBIN=/home/fit/Workspace/FiT/a.out
DATA=/mnt/data

## call bin
echo $HPARAMS | tr ' ' '\n' | parallel  -j $HJOBS "$FITBIN -f $DATA/$INPUT -o $DATA -h {} -g -k 100 -m1 -m2 -m3" &
echo $GPARAMS | tr ' ' '\n' | parallel  -j $GJOBS "$FITBIN -f $DATA/$INPUT -o $DATA -h -g {} -k 100 -m1 -m2 -m3" &
wait
