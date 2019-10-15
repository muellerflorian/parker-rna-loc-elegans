# Module with a bunch of general purpose utility functions.

# Imports
import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt

def cmap_create_random():
    """
    Random lookup table for Matplotlib. For instance, to show segmentation results.
    From: https://gist.github.com/jgomezdans/402500
    """

    vals = np.linspace(0,1,256)
    np.random.shuffle(vals)
    test = plt.cm.jet(vals)
    test[0,:] = [0.2, 0.2, 0.2, 0.2]   # Set to gray
    cmap = plt.cm.colors.ListedColormap(test)
    return cmap


def create_folder(folder_new):
    """
    Creates a folder if provided folder does not exist already

    Args:
        folder_new (string): name of folder to be created

    """

    if not os.path.isdir(folder_new):
        os.makedirs(folder_new)


# JSON encoder for numpy
# From https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def log_message(msg, callback_fun=None):
    """
    Output log, either terminal or some callback taking a string as input'''
    """
    if callback_fun:
        callback_fun(msg)
    else:
        print(msg)


def log_progress(iteration, total, callback_fun=None):
    """
    Progress log, either terminal or some callback taking a percantage as input''' 
    """
    if callback_fun:
        completed = (iteration / float(total))
        callback_fun(completed)
    else:
        
        print_progress(iteration, total)


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
    bar = '*' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s\r' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def file_name_decompose(file_name):
    """
    Decompose file name
    """

    drive, path_and_file = os.path.splitdrive(file_name)
    path, file = os.path.split(path_and_file)
    file_base, ext = os.path.splitext(file)

    return drive, path, file_base, ext


def create_mask_sphere(r):
    """
    Create a sphere with radius r.
    Returned array has size 2r+1 in all directions

    Args:
        r (integer): radius of sphere.

    Returns:
        mask (3D numpy array): binary mask of the sphere.
    """

    struct = np.zeros((2 * r + 1, 2 * r + 1, 2 * r + 1))
    x, y, z = np.indices((2 * r + 1, 2 * r + 1, 2 * r + 1))
    mask = (x - r) ** 2 + (y - r) ** 2 + (z - r) ** 2 <= r ** 2
    struct[mask] = 1
    return struct.astype('uint8')