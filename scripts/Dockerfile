FROM ubuntu:14.04

RUN apt-get update && apt-get install -y python2.7 python-pip python2.7-dev vim
RUN apt-get install -y --quiet wget 

RUN wget https://micans.org/mcl/src/mcl-12-068.tar.gz
RUN tar xzf mcl-12-068.tar.gz
WORKDIR "/mcl-12-068"
RUN ./configure
RUN make install
WORKDIR "/"

RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:webupd8team/java && apt-get update \
	&& echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | sudo /usr/bin/debconf-set-selections \
	&& apt-get install git curl -y oracle-java8-installer 

RUN apt-get install -y git
RUN git clone https://github.com/parthatalukdar/junto.git
ENV JUNTO_DIR="/junto"
ENV PATH="$PATH:$JUNTO_DIR/bin"

WORKDIR "/junto"
RUN bin/build update compile
RUN bin/build assembly
WORKDIR "/"

RUN pip install biogrouper
