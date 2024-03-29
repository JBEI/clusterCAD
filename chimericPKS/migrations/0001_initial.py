# Generated by Django 2.2.3 on 2020-03-09 21:34

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('pks', '0006_auto_20200309_2134'),
    ]

    operations = [
        migrations.CreateModel(
            name='ChimericPKS',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='Chimeric PKS', max_length=2000)),
                ('length', models.PositiveIntegerField()),
            ],
        ),
        migrations.CreateModel(
            name='ChimericSubunit',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='Chimeric Subunit', max_length=2000)),
                ('order', models.PositiveIntegerField(default='None')),
                ('length', models.PositiveIntegerField(default=0)),
                ('chimericPKS', models.ForeignKey(default=None, null=True, on_delete=django.db.models.deletion.CASCADE, to='chimericPKS.ChimericPKS')),
            ],
        ),
        migrations.CreateModel(
            name='ChimericFragment',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('order', models.PositiveIntegerField()),
                ('start', models.PositiveIntegerField()),
                ('stop', models.PositiveIntegerField()),
                ('front_sequences', models.TextField(default='unknown')),
                ('back_sequences', models.TextField(default='unknown')),
                ('length', models.PositiveIntegerField()),
                ('chimericSubunit', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='chimericPKS.ChimericSubunit')),
                ('naturalSubunit', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, to='pks.Subunit')),
            ],
        ),
    ]
