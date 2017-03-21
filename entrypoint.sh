#!/bin/bash

python3 manage.py migrate --noinput
python3 manage.py collectstatic --noinput
/usr/local/bin/gunicorn retrosyn.wsgi:application -w 2 -b :8000 --reload
