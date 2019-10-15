 # -*- coding: utf-8 -*-

"""
Introduction
============

Module containing different functions to work with FQ result files.

Usage
=====


"""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import math
import sys
import json
import numpy as np
import skimage
from nested_lookup import nested_lookup # pip install nested-lookup
from scipy import ndimage


# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------

# Info about the module
__author__ = "Florian MUELLER"
__email__ = "muellerf.research@gmail.com"


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

# Encoder class that allows to have numpy arrays in dictionary
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def colorbar(mappable):
    """
    Function to place colorbars next to images and guarantee that they have the
    same size.

    From: https://joseph-long.com/writing/colorbars/
    More info: https://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html#colorbar-whose-height-or-width-in-sync-with-the-master-axes
    """

    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)

    https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
    more info: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s\r' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def read_FQ_matlab(file_open):
    """ Opens FISH-quant result files generated with Matlab (tab-delimited text file).

    Args:
        file_open (string): string containing the full file name.

    Returns:
        dictionary containing outlines of cells, and if present the detected spots.
    """

    # Open file
    with open(file_open, "r") as fh:
         data = fh.readlines()

    # Strip white space characters
    data = [x.strip() for x in data]

     # Loop over read-in data
    fq_dict = {'cells':{},'file_names':{},'settings':{}}
    iLine = 0

    while iLine < len(data):

        line = data[iLine]


        # READ FILE NAMES
        if 'IMG_Raw' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'smFISH':img_name[1]})

        if 'IMG_Filtered' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'smFISH_filt':img_name[1]})

        if 'IMG_DAPI' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'DAPI':img_name[1]})

        if 'FILE_settings' in line:
            img_name = line.split('\t')
            if len(img_name) == 2:
                fq_dict['file_names'].update({'settings':img_name[1]})

        # READ IMAGE PARAMETERS
        if 'PARAMETERS' in line:
            iLine += 2
            par_microscope = data[iLine].split('\t')
            fq_dict['settings'].update({'microscope':{'pix_xy':float(par_microscope[0]),
                                                      'pix_z':float(par_microscope[1]),
                                                      'RI':float(par_microscope[2]),
                                                      'EX':float(par_microscope[3]),
                                                      'EM':float(par_microscope[4]),
                                                      'NA':float(par_microscope[5]),
                                                      'type':par_microscope[6]}})

        # New cell
        if 'CELL_START' in line:

            # Get name of cell
            cell_id = line.split('\t')[1]

            ### POSITION OF CELL

            # Read X-POS
            iLine += 1
            pos_list = (data[iLine].replace('X_POS\t','')).split('\t')
            x_pos    = [int(s) for s in pos_list]

            # Read Y-POS
            iLine += 1
            pos_list = (data[iLine].replace('Y_POS\t','')).split('\t')
            y_pos    = [int(s) for s in pos_list]

            # Read Z-POS
            iLine += 1
            pos_list = (data[iLine].replace('Z_POS\t','')).split('\t')
            if len(pos_list) > 1:
                z_pos = [int(s) for s in pos_list]
            else:
                z_pos = ['']

            fq_dict['cells'].update({cell_id:{'cell_pos':{'x': x_pos,'y': y_pos,'z': z_pos}}})


        # New nucleus
        if 'Nucleus_START' in line:

            # Get name of cell
            nuc_id = line.split('\t')[1]

            ### POSITION OF CELL

            # Read X-POS
            iLine += 1
            pos_list = (data[iLine].replace('X_POS\t','')).split('\t')
            x_pos    = [int(s) for s in pos_list]

            # Read Y-POS
            iLine += 1
            pos_list = (data[iLine].replace('Y_POS\t','')).split('\t')
            y_pos    = [int(s) for s in pos_list]

            # Read Z-POS
            iLine += 1
            pos_list = (data[iLine].replace('Z_POS\t','')).split('\t')
            if len(pos_list) > 1:
                z_pos = [int(s) for s in pos_list]
            else:
                z_pos = ['']

            fq_dict['cells'][cell_id].update({nuc_id:{'nuc_pos':{'x': x_pos,'y': y_pos,'z': z_pos}}})



        # Position of detected RNAS
        if 'SPOTS_START' in line:
            iLine += 2 # Move over header

            RNA_prop = []
            while not('SPOTS_END' in data[iLine]):
                RNA_prop.append([float(s) for s in data[iLine].split('\t')])
                iLine += 1

            # Assign to dictionary
            fq_dict['cells'][cell_id].update({'spots': np.array(RNA_prop)})


        # Up date line counter
        iLine += 1

    return fq_dict



def get_rna(fq_dict):
    """
    Obtain a numpy array with all detected spots in the image. Detection results
    are saved in a dictionary (see read_FQ_results_matlab for more details).
    """
    RNAall = nested_lookup('spots', fq_dict)  # returns list of numpy arrays

    for idx,val in enumerate(RNAall):
        if idx == 0:
            spots_all = np.copy(val)
        else:
            spots_all = np.append(spots_all,val,axis=0)

    return spots_all


def calc_expression_density_plot(fq_dict,img_size,outline_int = 'max',flag_plot = False):
    """ Calculate expression density image.
    RNA detection results are used to calculate a 2D image where each cell
    is displayed with pixel values corresponding to the number of RNAs in the cell.

    Args:
        imageprop ('dict'): dictionary containing information about outlines of cells
        and nuclei as well as (if present) positions of RNA molecules

        img_size (tuple): specifying the size of the image.

        outline_int (string) specifying how pixel values of cell outlines in the
        density plot.'max' means that the maximum number of RNAs per cell is used.
        '*nt' with int being the integer value that should be used.

        flag_plot ('bool'): flag to indicate if results should be plotted.

    Returns:
        2D numpy arrays (i) image with outlines, (ii) image with
        expression density, (iii) image wiht expression density and outlines.

    """

    img_density = np.zeros(img_size, dtype=np.uint16)
    img_outline = np.zeros(img_size, dtype=np.uint8)

    # Generate image of outline and density
    iCell = 1
    print_progress(iCell, len(fq_dict['cells']))

    for key, value in fq_dict['cells'].items():
        print_progress(iCell, len(fq_dict['cells']))

        cell_pos = []
        cell_pos.append(value['cell_pos']['x'])
        cell_pos.append(value['cell_pos']['y'])
        cell_pos = np.array(cell_pos)

        # How many RNAs
        if 'spots' in value.keys():
            Nrna = value['spots'].shape[0]
        else:
            Nrna = 0

        # Create contour image
        [rr_cont, cc_cont] = skimage.drawSK.polygon_perimeter(cell_pos[1,:], cell_pos[0,:], shape=img_outline.shape, clip=True)
        img_outline[rr_cont, cc_cont] = 1

        # Generate coordinates of pixels within polygon.
        rr, cc = skimage.drawSK.polygon(cell_pos[1,:], cell_pos[0,:])
        img_density[rr, cc] = Nrna

        # Update cell counter
        iCell += 1


    ## Add outline mask to density plot

    # Decide which intensity the outline should have
    if outline_int == 'max':
        line_int = np.amax(img_density)
    else:
        line_int = float(outline_int)

    img_density_outline  = np.copy(img_density)
    img_density_outline[img_outline==1] = line_int


    # Plot results
    if flag_plot:
        fig, (ax1, ax2,ax3) = plt.subplots(3,1,num='density_plt')
        img1 = ax1.imshow(img_density,cmap="hot")
        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)
        colorbar(img1)

        img2 = ax2.imshow(img_density_outline,cmap="hot")
        ax2.get_xaxis().set_visible(False)
        ax2.get_yaxis().set_visible(False)
        colorbar(img2)

        ax3.imshow(img_outline,cmap="hot")
        ax3.get_xaxis().set_visible(False)
        ax3.get_yaxis().set_visible(False)

        plt.tight_layout(h_pad=1)
        plt.draw()

    # Return results
    return img_density,img_density_outline,img_outline
