#!/bin/bash

python3 manage.py migrate --noinput
python3 manage.py collectstatic --noinput
/usr/local/bin/gunicorn clusterCAD.wsgi:application -w 5 -b :8000 --reload
