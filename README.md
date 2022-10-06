---------------------
FiT anomaly detection
---------------------

Description
===========

Given a link stream in input, aggregate the stream into different graph sizes and compute features for each interaction.

Installation
============

Requirements
------------

You need a g++ compilator compatible with C++14, cmake (>=3.10), networkit and catch2 ( https://github.com/catchorg/Catch2 )

NetworKit installation
----------------------

First you need to install networkit as a library :

```
mkdir </PATH/TO/NETWORKIT/LIB/> # choose a directory where you will install networkit library
git clone --recursive https://github.com/networkit/networkit.git
cd networkit
mkdir build && cd build

cmake -DCMAKE_INSTALL_PREFIX:PATH=</PATH/TO/NETWORKIT/LIB/> ..
make -j5 && make install
```
An error might occur where `ttmath` headers are not copied in the library. You need to check if `</PATH/TO/NETWORKIT/LIB/>/include/ttmath/` is empty.
If it is empty, copy the content of `networkit/extlibs/ttmath/ttmath/` into `</PATH/TO/NETWORKIT/LIB/>/include/ttmath/`

Then 

```
export LD_LIBRARY_PATH=</PATH/TO/NETWORKIT/LIB/>/lib
```

Project compilation
-------------------

When NetworKit is installed, you can compile the project:

```
mkdir build && cd build
cmake ../
make 
```

You can then run the unit tests to ensure that everything works as intended:

```
./tests
```

Usage
-----

The parameter order is important.
The parameters noted as "required" are mandatory, the "optional"  are not.
: 
```
    ./fit.bin  [required] -f <Path To CSV> -o <Path To Output Folder> -h <List Of H Sizes> -g <List Of G Sizes> -k <Bound On The Degree> [optional] -m1 -m2 -m3 -p -b
```

Example (for 1 run with 1 history graph H of size 1024) :
```
    ./fit.bin  -f uci/uci.010/added_links.txt -o uci_output/ -h 1024 -g -k 10  -m1 -m2 -p
```

Optional parameters [see below for list of metrics toggled by m1, m2 and m3): 

-m1 : toggles basic metrics

-m2 : toggles local metrics (in O(k), O(k²))

-m3 : toggles more intensive metrics

-p : when enabled, compute projection graph, and output metrics on them

-b : enable when the input graph is bipartite



Using Docker (DEPRECATED - TODO: update) 
------------

Build the docker image using 

  docker build -t getting-started .

To use the docker, you then have to create a data/ folder inside your current directory, put your <file> in data/, then run using 

  docker run -v $PWD/data:/mnt/data -it getting-started <file> <njobs>

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
    * link component number of links
    * number of components in graph
    * size of largest component
    * eccentricity of u in G
    * eccentricity of v in G
    * eccentricity of u in G-
    * eccentricity of v in G-
    * total distance change u
    * total distance change v
    * max distance change u
    * max distance change v
    * number distance change u
    * number distance change v
    * sum difference of distances to u and v in G
    * sum difference of distances to u and v in G-
    * closeness of u in G
    * closeness of v in G
    * closeness of u in G-
    * closeness of v in G-
    * distance u v in G-
    * degeneracy G
    * degeneracy G-
    * sum core scores G
    * sum core scores G-
    * core number u in G
    * core number v in G
    * core number u in G-
    * core number v in G-
    * sum core variations
    * number of components in G-
    * size largest component G-
    * pagerank of u in G
    * pagerank of v in G
    * pagerank of u in G-
    * pagerank of v in G-
    * pagerank max in G
    * pagerank max in G-
    * pagerank variation
    * number of subsets in G-
    * max subset size in G-
    * u subset size in G-
    * v subset size in G-
    * u v in same community in G-
    * number of subsets in G
    * max subset size in G
    * u subset size in G
    * v subset size in G
    * u v in same community in G
    * number of nodes changing partition
