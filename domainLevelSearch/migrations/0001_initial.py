# Generated by Django 2.2.4 on 2022-09-14 19:48

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('pks', '0002_auto_20220914_1948'),
    ]

    operations = [
        migrations.CreateModel(
            name='DomainChar',
            fields=[
                ('id', models.BigAutoField(primary_key=True, serialize=False)),
                ('domainString', models.CharField(db_index=True, max_length=1000, unique=True)),
            ],
        ),
        migrations.CreateModel(
            name='ClusterString',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('archString', models.CharField(db_index=True, max_length=1000)),
                ('cluster', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='pks.Cluster')),
            ],
        ),
    ]
