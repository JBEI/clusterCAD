#!/bin/bash

# this script launches this service on an aws container, after manually adding the servers ssh key to bitbucket
# make sure to enable port 80 in Security groups for this container as well
# also make sure you are using DEBUG = False, and settings ALLOWED_HOSTS to the correct server name
# for production

sudo yum update -y
sudo yum install -y docker git
sudo pip install docker-compose
sudo service docker start

# create ssh key and add to bitbucket before running this
git clone ssh://git@repo.jbei.org:7999/pks/pksretrosyn.git
cd pksretrosyn
sudo /usr/local/bin/docker-compose build
sudo /usr/local/bin/docker-compose up -d

# example to forward ports back for devel
ssh -L 80:localhost:80 8888:localhost:8888 ec2-user@<remote hostname>
