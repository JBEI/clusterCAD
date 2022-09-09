from django.db import models
import pks.models

class DomainChar(models.Model):
    # This class stores a unique integer for each type of domain in
    # ClusterCAD, each represented by a unique text string

    id = models.BigAutoField(primary_key=True)
    domainString = models.CharField(max_length=1000, db_index=True, blank=False)

    def char(self):
        # returns the unique unicode char for this domain
        return chr(id)

class ClusterString(models.Model):
    # This class stores a string representing each PKS architecture
    # with one domain per string

    cluster = models.OneToOneField(pks.models.Cluster, on_delete=models.CASCADE)
    archString = models.CharField(max_length=1000, db_index=True, blank=False)
