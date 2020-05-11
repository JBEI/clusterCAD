import sys
from django.db import models
from model_utils.managers import InheritanceManager
from rdkit import Chem as chem
from rdkit.Chem import AllChem, rdFMCS
import os
import json
from collections import OrderedDict
from compounddb.models import Compound
from django.db.models.signals import pre_save, pre_delete
from django.dispatch import receiver
from django.core.validators import validate_comma_separated_integer_list
from pdb import set_trace as st
import numpy as np
# from pks.models import Subunit


#Question: bottom up, lower structures contain reference to upper layers, then subunit should have reference to chimeric fragments?

class ChimericPKS(models.Model):

    # A Chimeric PKS must has certain structures to be a valid pks, or the chemistry won't succeed.
    # This is the class of the chimericPKS, which has chimericSubunits as the sub-structure.
    '''
    Multiple chimeric subunits, each representing a single peptide/protein/gene with one or more PKS modules
    Each subunit represents cut together sections from parent subunits, of natural PKSs
    Point mutations can be made, as annotations to the sequences of the subunits
    Coordinates at the amino acid level, with ability to export either original DNA, or create codon optimized DNA
    Export to J5/DIVA or GenBank format for building PKSs
    '''
    '''Properties:
    name: name of this chimeric PKS (default chimeric PKS)
    length: the length of the PKS structure, in terms of number of chimeric subunits
    '''


    name = models.CharField(max_length=2000, default='Chimeric PKS')
    length = models.PositiveIntegerField()

    # Returns a list of chimeric subunits 
    def subunits(self):
        return ChimericSubunit.objects.filter(chimericPKS=self).order_by('order')

    # Upon calling, sorts the subunits so that order is incremental by 1, starting at 0
    def reorder(self):
        counter = 0
        for subunit in self.subunits():
            subunit.order = counter
            counter+=1
            subunit.save()



    # def insert(self, positions, fragments):
    #     #Support add last
    #     length = len(self.subunits())
    #     last=False
    #     if positions[-1] == length+1:
    #         last=True
    #         positions = positions[:-1]
            
    #     total_length = length
    #     old_fragment_update=[]
    #     position_index = 0
    #     previous_index = 0
    #     added_term = 0
    #     previous = np.arange(length)
    #     new_index = []
    #     while position_index + previous_index < len(positions) + length:
    #         if position_index == len(positions) or previous[previous_index] + added_term < positions[position_index] + min(0, (added_term-1)*-1)*-1:
    #             old_fragment_update.append(previous[previous_index] + added_term)
    #             previous_index += 1
    #         else:
    #             # old_fragment_update.append(positions[position_index] + min(0, (added_term-1)*-1)*-1)
    #             new_index.append(positions[position_index] + min(0, (added_term-1)*-1)*-1)
    #             position_index += 1
    #             added_term += 1
    #     if last:
    #         last_value = position_index + previous_index 
    #         # old_fragment_update.append(last_value)
    #         new_index.append(last_value)

    #     st()
    #     index=0
    #     for frag in self.subunits():
    #         frag.order = old_fragment_update[index]
    #         index+=1
    #         frag.save()

    #     index=0
    #     for frag in fragments:
    #         frag.order = new_index[index]
    #         index+=1
    #         frag.save()

    # Given subunits, which is a list of names, 
    # remove all instances within that subunit with that name
    def delete(self, subunits):
        # Subunits are a list of names
        for subunit in subunits:
            ChimericSubunit.objects.filter(name=subunit.name, chimericPKS=self).delete()
        self.reorder()

    # Remove a subunit specified by the order.
    def delete_on_order(self, order):
        # Subunits are a list of names
        ChimericSubunit.objects.filter(order=order, chimericPKS=self).delete()
        self.reorder()

    # def exportJ5(self):
    #   return

    # def __str__(self):
    #     return "%s gene cluster" % self.description

class ChimericSubunit(models.Model):

    # This is the class for chimeric subunits, which is affilated to chimericPKS, and has substructure chimericFragments
    '''
    list of fragments, fragments themselves should be separate objects (each fragment is a piece of information from a subunit)
    Properties:
        name: name of this chimeric Subunit (default chimeric Subunit)
        chimericPKS: the parent chimericPKS that contains this subunit
        order: the order in which this subunit is stored in the parent chimericPKS structure
        length: the number of chimeric fragments this contains
    '''
    name = models.CharField(max_length=2000, default='Chimeric Subunit')
    chimericPKS = models.ForeignKey(ChimericPKS, default=None, on_delete=models.CASCADE, null=True)
    order = models.PositiveIntegerField(default='None')
    length = models.PositiveIntegerField(default=0)
    # visible = models.BooleanField(default=True)

    # This is code for printing the subunit structure, currently it is of the form
    # NAME_order: **position within chimericPKS**
    def __repr__(self):
        return self.name+'_order:'+str(self.order)

    # Gets the fragments in order
    def get_order(self):
        return ChimericFragment.objects.filter(chimericSubunit=self).order_by('order')

    # Reordering
    def reorder(self):
        counter = 0
        for fragment in self.fragments():
            fragment.order = counter
            counter+=1

    # This function is no longer used, as insertion is handled upon creation
    # def insert_some(self, positions, fragments):
    #     #Support add last
    #     length = self.length
    #     last=False
    #     if positions[-1] == length+1:
    #         last=True
    #         positions = positions[:-1]
            
    #     total_length = length
    #     old_fragment_update=[]
    #     position_index = 0
    #     previous_index = 0
    #     added_term = 0
    #     previous = np.arange(length)
    #     new_index = []
    #     while position_index + previous_index < len(positions) + length:
    #         if position_index == len(positions) or previous[previous_index] + added_term < positions[position_index] + min(0, (added_term-1)*-1)*-1:
    #             old_fragment_update.append(previous[previous_index] + added_term)
    #             previous_index += 1
    #         else:
    #             # old_fragment_update.append(positions[position_index] + min(0, (added_term-1)*-1)*-1)
    #             new_index.append(positions[position_index] + min(0, (added_term-1)*-1)*-1)
    #             position_index += 1
    #             added_term += 1
    #     if last:
    #         last_value = position_index + previous_index 
    #         # old_fragment_update.append(last_value)
    #         new_index.append(last_value)
    #     index=0
    #     for frag in self.get_order():
    #         frag.order = old_fragment_update[index]
    #         index+=1

    #     index=0
    #     for frag in fragments:
    #         frag.order = new_index[index]
    #         index+=1

    def delete(self, fragments):
        # Subunits are a list of names
        for fragment in fragments:
            ChimericFragment.objects.filter(name=fragment, chimericSubunit=self).delete()
        self.reorder()

    # Remove a fragment specified by the order.
    def delete_on_order(self, order):
        # Subunits are a list of names
        ChimericFragment.objects.filter(order=order, chimericSubunit=self).delete()
        self.reorder()


@receiver(pre_save, sender=ChimericSubunit)
# This function gets called everytime we call .save() for structures, 
# this will sort the items in correct order automatically
def setSubunitOrder(sender, instance, **kwargs):
    # sets the subunit order
    # If order has [0, 1, 2, 3] and input 4, then get [0, 1, 2, 3] + [4]
    # If order has [0, 1, 2, 3] and input 3, then get [0, 1, 2] + [3] + [4]
    if not isinstance(instance.order, int):
        try:
            subunitCount = sender.objects.filter(chimericPKS=instance.chimericPKS).count()
            instance.order = subunitCount
        except:
            instance.order = 999999
            assert 'Something wrong'
            return
    else:
        for subunit in list(sender.objects.filter(chimericPKS=instance.chimericPKS)):
            if subunit.order == instance.order:
                subunit.order = subunit.order + 1
                subunit.save()
                break
            


    #add fragment to set order
    
    # def fragments(self):
    #     return Chimeric.objects.filter(cluster=self).order_by('order')

    # def annotate(self):
    #   return 

class ChimericFragment(models.Model):
    # I haven't coded much for this class, this class essentially inherits from real molecule fragments
    chimericSubunit = models.ForeignKey(ChimericSubunit, on_delete=models.CASCADE)
    naturalSubunit = models.ForeignKey('pks.Subunit', on_delete=models.DO_NOTHING)

    order = models.PositiveIntegerField()
    #Amino acid coordinates
    start = models.PositiveIntegerField()
    stop = models.PositiveIntegerField()
    front_sequences = models.TextField(default='unknown')
    back_sequences = models.TextField(default='unknown')
    length = models.PositiveIntegerField()


