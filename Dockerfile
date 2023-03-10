FROM ubuntu:22.04

Label maintainer="rdemko2332@gmail.com"

WORKDIR /usr/bin/

RUN apt-get -qq update --fix-missing \
  && apt-get install -y wget perl libgomp1 git ant build-essential unzip default-jre python3 cpanminus bioperl libaio1 emacs libjson-perl libmodule-install-rdf-perl libxml-parser-perl openjdk-8-jdk libdate-manip-perl libtext-csv-perl libstatistics-descriptive-perl libtree-dagnode-perl libxml-simple-perl && apt-get clean \
  && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*

RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.4/diamond-linux64.tar.gz
RUN tar xzf diamond-linux64.tar.gz

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

ADD /bin/*.pl /usr/bin/

# Making all tools executable
RUN chmod +x *

WORKDIR /work