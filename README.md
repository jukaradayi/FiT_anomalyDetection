---------------------
FiT anomaly detection
---------------------

Description
-----------

Given a link stream in input, aggregate the stream into different graph sizes and compute features for each interaction.

Usage
-----

The parameter order is important.
The parameters noted as "required" are mandatory, the "optional"  are not.
: 
```
    ./a.out  [required] -f <Path To CSV> -o <Path To Output Folder> -h <List Of H Sizes> -g <List Of G Sizes> -k <Bound On The Degree> [optional] -m1 -m2 -m3 -p -b
```

Example :
```
    ./a.out  -f uci/uci.010/added_links.txt -o uci_output/ -h 1024 2048 4096 -g 60 3600 7200  -k 10  -m1 -m2 -p
```

Optional parameters: 

-m1 : toggles basic metrics

-m2 : toggles local metrics (in O(k), O(k²))

-m3 : toggles more intensive metrics

-p : when enabled, compute projection graph, and output metrics on them

-b : enable when the input graph is bipartite

Requirement and Installation
----------------------------

The only requirement is to build networkit (any version, should work with the latest) as a library.

- A modern C++ compiler (g++ >= 4.8) compatible with c++11
- OpenMP (usually ships with the compiler)

To build networkit : 

```
mkdir build
cd build
cmake ..
make -jX install
```

Then to build the project: 

```
g++ -L<PATH_TO_NETWORKIT>/lib -I<PATH_TO_NETWORKIT>/include -lnetworkit -std=c++11 -fopenmp main.cpp metrics.cpp G_graph.cpp history_graph.cpp H_graph.cpp PLM.cpp

```
