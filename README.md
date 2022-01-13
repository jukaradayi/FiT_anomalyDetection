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

Optional parameters [see below for list of metrics toggled by m1, m2 and m3): 

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

List of Metrics
---------------
(metrics including top need "-p" to be outputed, metrics including bot need "-p -b" to be outputed)
* m1: 
    * number of nodes
    * number of links
    * number of nodes of degree 1
    * number of nodes of degree 2
    * degree of u
    * degree of v
    * weighted degree of u
    * weighted degree of v
    * max degree
    * max weighted degree
    * degree absolute difference
    * weighted degree absolute difference
    * weight of link u-v
    * total weight of the graph
    * average weight
    * average degree
    * average weighted degree
    * density
    * number of nodes in top graph (if -p enabled)
    * number of nodes in bot graph (if -p -b)
    * degree of u in top graph (if -p)
    * degree of v in bot graph (if -b -p)
    * top weighted degree (if -p)
    * bot weighted degree (if -p -b)
    * top max weighted degree (if -p)
    * bot max weighted degree (if -p -b)

* m2:
    * Jaccard bipartite u
    * Jaccard bipartite v
    * ratio of neighbors of u in N(v)
    * ratio of neighbors of v in N(u)
    * number of common neighbors expected from degree
    * neighborhood overlap
    * egonet Nu number of links
    * egonet Nu number of nodes
    * egonet Nv number of links
    * egonet Nv number of nodes
    * egonet Nu u Nv number of links
    * egonet Nu u Nv number of nodes
    * egonet Nu u NNu number of links
    * egonet Nu u NNu number of nodes
    * egonet Nv u NNv number of links
    * egonet Nv u NNv  number of nodes
    * egonet (Nu u NNv) u (NNu u Nv) number of links
    * egonet (Nu u NNv) u (NNu u Nv) number of nodes
    * egonet (Nu u NNv) n (NNu u Nv) number of links
    * egonet (Nu u NNv) n (NNu u Nv) number of nodes
    * egonet Nu u Nv maxsize
    * egonet Nv u NNv maxsize
    * egonet Nu u NNu maxsize
    * egonet (NNu u Nv) u (Nu u NNv) maxsize
    * egonet (NNu u Nv) n (Nu u NNv) maxsize
    * addamic adar
    * adamic adar bipartite u 
    * adamic adar bipartite v
    * top clustering u
    * bot clustering v
* m3 : 
    * link component number of nodes
    *  link component number of links



    
