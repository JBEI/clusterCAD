# -*- coding: utf-8 -*-
# Generated by Django 1.10.5 on 2017-03-14 23:36
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pks', '0004_subunit_genbankaccession'),
    ]

    operations = [
        migrations.AlterField(
            model_name='subunit',
            name='genbankAccession',
            field=models.CharField(max_length=2000),
        ),
    ]
