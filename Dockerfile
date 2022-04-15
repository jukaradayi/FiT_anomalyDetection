FROM ubuntu:bionic

USER root

RUN apt-get update -y; \
apt-get install -y build-essential git cmake autoconf libtool pkg-config cmake vim wget parallel;

RUN useradd -ms /bin/bash fit;

#USER fit

RUN mkdir -p  /home/fit/Workspace ; \
mkdir -p /home/fit/Workspace/FiT; \
cd /home/fit/Workspace && git clone -b 9.1 --recursive https://github.com/networkit/networkit.git && mkdir networkit_lib ; \
cd /home/fit/Workspace/networkit && mkdir build && cd build ; \ 
cmake -DCMAKE_INSTALL_PREFIX:PATH=/home/fit/Workspace/networkit_lib .. && make -j 5 && make install ; 
#cp /home/fit/Workspace/networkit/extlibs/ttmath/ttmath/*hpp /home/fit/Workspace/networkit_lib/include/ttmath/; # for some reason ttmath headers are not copied in install dir


ADD ./*pp /home/fit/Workspace/FiT/
ADD ./run.sh /home/fit/Workspace/FiT/


RUN cd /home/fit/Workspace/FiT/ ;\
chmod +x run.sh; \
ls; 

RUN export LD_LIBRARY_PATH=/home/fit/Workspace/networkit_lib/; \
cd /home/fit/Workspace/FiT/ && g++ main.cpp metrics.cpp G_graph.cpp history_graph.cpp H_graph.cpp PLM.cpp -L/home/fit/Workspace/networkit_lib/lib -I/home/fit/Workspace/networkit_lib/include -lnetworkit -fopenmp -std=c++14 -O3;

ENTRYPOINT [ "/home/fit/Workspace/FiT/run.sh"]
