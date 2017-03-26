# -*- coding: utf-8 -*-
# Generated by Django 1.10.5 on 2017-03-26 22:15
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('compounddb', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Cluster',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('genbankAccession', models.CharField(max_length=2000, unique=True)),
                ('mibigAccession', models.CharField(max_length=2000, unique=True)),
                ('description', models.TextField()),
                ('sequence', models.TextField()),
                ('knownProductSource', models.TextField()),
                ('knownProduct', models.ForeignKey(blank=True, default=None, null=True, on_delete=django.db.models.deletion.SET_NULL, to='compounddb.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='Domain',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('start', models.PositiveIntegerField()),
                ('stop', models.PositiveIntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='Module',
            fields=[
                ('order', models.AutoField(primary_key=True, serialize=False)),
                ('loading', models.BooleanField()),
                ('terminal', models.BooleanField()),
                ('product', models.ForeignKey(blank=True, default=None, null=True, on_delete=django.db.models.deletion.SET_NULL, to='compounddb.Compound')),
            ],
        ),
        migrations.CreateModel(
            name='Standalone',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('order', models.PositiveSmallIntegerField()),
                ('name', models.CharField(max_length=2000)),
                ('start', models.PositiveIntegerField()),
                ('stop', models.PositiveIntegerField()),
                ('sequence', models.TextField()),
                ('cluster', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='pks.Cluster')),
            ],
        ),
        migrations.CreateModel(
            name='Subunit',
            fields=[
                ('order', models.AutoField(primary_key=True, serialize=False)),
                ('genbankAccession', models.CharField(max_length=2000)),
                ('name', models.CharField(max_length=2000)),
                ('start', models.PositiveIntegerField()),
                ('stop', models.PositiveIntegerField()),
                ('sequence', models.TextField()),
                ('cluster', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='pks.Cluster')),
            ],
        ),
        migrations.CreateModel(
            name='ACP',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='AT',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
                ('substrate', models.CharField(choices=[('mal', 'mal'), ('mmal', 'mmal'), ('mxmal', 'mxmal'), ('emal', 'emal'), ('cemal', 'cemal'), ('Acetyl-CoA', 'Acetyl-CoA'), ('prop', 'prop'), ('isobut', 'isobut'), ('2metbut', '2metbut'), ('CHC-CoA', 'CHC-CoA'), ('trans-1,2-CPDA', 'trans-1,2-CPDA'), ('N/A', 'N/A')], default='mal', max_length=20)),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='cMT',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
                ('active', models.BooleanField(default=True)),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='DH',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
                ('active', models.BooleanField(default=True)),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='ER',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
                ('active', models.BooleanField(default=True)),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='KR',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
                ('active', models.BooleanField()),
                ('type', models.CharField(choices=[('A1', 'A1'), ('A2', 'A2'), ('B1', 'B1'), ('B2', 'B2'), ('C1', 'C1'), ('C2', 'C2'), ('U', 'U')], default=None, max_length=2)),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='KS',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='oMT',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
                ('active', models.BooleanField(default=True)),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='PCP',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
            ],
            bases=('pks.domain',),
        ),
        migrations.CreateModel(
            name='TE',
            fields=[
                ('domain_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='pks.Domain')),
                ('cyclic', models.BooleanField()),
            ],
            bases=('pks.domain',),
        ),
        migrations.AddField(
            model_name='module',
            name='subunit',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='pks.Subunit'),
        ),
        migrations.AddField(
            model_name='domain',
            name='module',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='pks.Module'),
        ),
    ]
