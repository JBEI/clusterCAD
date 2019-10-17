FROM debian:10

# Maintained by Tyler W. H. Backman
LABEL maintainer="tbackman@lbl.gov"

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
	    python3-tk \
        postgresql \
        postgresql-client \
	    postgresql-server-dev-all \
	    postgresql-server-dev-11 \
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
        libboost-all-dev \
        libboost-system1.67-dev \
        libboost-thread1.67-dev \
        libboost-serialization1.67-dev \
        libboost-python1.67-dev \
        libboost-regex1.67-dev \
        libboost-iostreams1.67-dev \
        libxrender1 \
        libxext6 \ 
        libeigen3-dev \
        openjdk-11-jdk \
        openjdk-11-jre \
        sudo \
        ncbi-blast+ \
	    libigraph0v5 \
	    libigraph0-dev \
        zlib1g-dev \
        libcairo2-dev \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# install pip packages from requirements.txt
# note: requirements.txt can be updated with 'pip3 install pur; pur -r requirements.txt'
RUN pip3 install -U pip setuptools
RUN pip3 install -r requirements.txt


# compile and install rdkit from source code
ENV RDKIT_VERSION Release_2019_03_3
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
        -D Py_ENABLE_SHARED=1 \
        -D PYTHON_EXECUTABLE=/usr/bin/python3 \
        -D RDK_INSTALL_INTREE=ON \
        -D RDK_BUILD_INCHI_SUPPORT=ON \
        -D RDK_BUILD_AVALON_SUPPORT=ON \
        -D RDK_BUILD_PYTHON_WRAPPERS=ON \
        -D RDK_BUILD_CAIRO_SUPPORT=ON \
        -D RDK_BUILD_CPP_TESTS=ON \
	    -D RDK_BUILD_PGSQL=ON \
	    -D PostgreSQL_TYPE_INCLUDE_DIR=/usr/include/postgresql/11/server \
	    -D PostgreSQL_ROOT=/usr \
        .. && \
    make -j `grep -c ^processor /proc/cpuinfo` && \
    make install 
RUN chmod u+x /root/rdkit-$RDKIT_VERSION/build/Code/PgSQL/rdkit/pgsql_install.sh && \
    /root/rdkit-$RDKIT_VERSION/build/Code/PgSQL/rdkit/pgsql_install.sh

# create postgresql user and database
# and install RDKit cartridge
RUN echo "host all  all    0.0.0.0/0  md5" >> /etc/postgresql/11/main/pg_hba.conf
RUN echo "listen_addresses='*'" >> /etc/postgresql/11/main/postgresql.conf
RUN echo "shared_buffers = 1024MB" >> /etc/postgresql/11/main/postgresql.conf
RUN echo "work_mem = 128MB" >> /etc/postgresql/11/main/postgresql.conf
RUN echo "synchronous_commit = off" >> /etc/postgresql/11/main/postgresql.conf
USER postgres
WORKDIR /
RUN /etc/init.d/postgresql start &&\
    psql -c "CREATE USER clustercad WITH SUPERUSER PASSWORD 'clustercad';" &&\
    createdb -O clustercad clustercad && \
    psql -d clustercad -c "create extension if not exists rdkit; create schema rdk;"
VOLUME  ["/etc/postgresql", "/var/log/postgresql", "/var/lib/postgresql"]
EXPOSE 5432

USER root
