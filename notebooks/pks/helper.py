import numpy as np
from matplotlib import pyplot as plt
from rdkit.Chem import Draw



def plot_mol(mol):
    '''Helper function to plot molecule on new figure.
    '''
    plt.imshow(np.asarray(Draw.MolsToImage([mol])))

