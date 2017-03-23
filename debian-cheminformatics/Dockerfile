FROM debian:8 

MAINTAINER Tyler W. H. Backman "tbackman@lbl.gov"

USER root

ENV LANG C.UTF-8
WORKDIR /root
ADD ./requirements.txt /root/

# include Debian packages required to build pip packages
RUN apt-get update \
    && apt-get install -y \
        python3-pip \
        python3-all \
        python3-all-dev \
        postgresql \
        postgresql-client \
	postgresql-server-dev-all \
	postgresql-server-dev-9.4 \
        libpq-dev \
        vim \
        wget \
        unzip \
        build-essential \
        cmake \
        python3-dev \
        sqlite3 \
        libsqlite3-dev \
        libboost-dev \
        libboost-system-dev \
        libboost-thread-dev \
        libboost-serialization-dev \
        libboost-python-dev \
        libboost-regex-dev \
        libxrender1 \
        libxext6 \ 
        libeigen3-dev \
        openjdk-7-jdk \
        openjdk-7-jre \
        sudo \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# install pip packages from requirements.txt
RUN pip3 install -U pip setuptools
RUN pip3 install -r requirements.txt

# compile and install rdkit from source code
ENV RDKIT_VERSION Release_2016_09_3
ENV RDBASE /root/rdkit-$RDKIT_VERSION
ENV LD_LIBRARY_PATH $RDBASE/lib
ENV PYTHONPATH $PYTHONPATH:$RDBASE
RUN wget https://github.com/rdkit/rdkit/archive/$RDKIT_VERSION.tar.gz && \
    tar xfz $RDKIT_VERSION.tar.gz && \
    rm $RDKIT_VERSION.tar.gz && \
    cd /root/rdkit-$RDKIT_VERSION/External/INCHI-API && \
    ./download-inchi.sh && \
    cd /root/rdkit-$RDKIT_VERSION && \
    mkdir build && \
    cd build && \
    cmake \
        -DRDK_BUILD_INCHI_SUPPORT=ON \
        -D PYTHON_LIBRARY=/usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4m.a \
        -D PYTHON_INCLUDE_DIR=/usr/include/python3.4/ \
        -D PYTHON_EXECUTABLE=/usr/bin/python3 \
	-D RDK_BUILD_PGSQL=ON \
	-D PostgreSQL_TYPE_INCLUDE_DIR=/usr/include/postgresql \
	-D PostgreSQL_ROOT=/usr \
        .. && \
    make -j `grep -c ^processor /proc/cpuinfo` && \
    make install 
RUN /root/rdkit-$RDKIT_VERSION/build/Code/PgSQL/rdkit/pgsql_install.sh

# create postgresql user and database
# and install RDKit cartridge
RUN echo "host all  all    0.0.0.0/0  md5" >> /etc/postgresql/9.4/main/pg_hba.conf
RUN echo "listen_addresses='*'" >> /etc/postgresql/9.4/main/postgresql.conf
RUN echo "shared_buffers = 1024MB" >> /etc/postgresql/9.4/main/postgresql.conf
RUN echo "work_mem = 128MB" >> /etc/postgresql/9.4/main/postgresql.conf
RUN echo "synchronous_commit = off" >> /etc/postgresql/9.4/main/postgresql.conf
USER postgres
WORKDIR /
RUN /etc/init.d/postgresql start &&\
    psql --command "CREATE USER clustercad WITH SUPERUSER PASSWORD 'clustercad';" &&\
    createdb -O clustercad clustercad
VOLUME  ["/etc/postgresql", "/var/log/postgresql", "/var/lib/postgresql"]
EXPOSE 5432

USER root
