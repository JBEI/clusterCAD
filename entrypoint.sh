#!/bin/bash

# copy example settings.py into place only if it doesn't already exist
# this allows production to maintain it's own unique settings.py
if [ ! -f clusterCAD/settings.py ]; then
    echo "settings.py not found, copying in example file"
    cp clusterCAD/example_settings.py clusterCAD/settings.py
else
    echo "existing settings.py file found, delete this first if you want to use the example version"
fi

python3 manage.py migrate --noinput
python3 manage.py collectstatic --noinput
/usr/local/bin/gunicorn clusterCAD.wsgi:application -w 5 -b :8000 --reload
