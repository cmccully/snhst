import shutil
import os


def copy_if_not_exists(filename, path):
    mkdir(path)
    copied_image = os.path.join(path, os.path.basename(filename))
    if not os.path.exists(copied_image):
        shutil.copy(filename, path)
        print(filename + ' --> ' + path)
    return copied_image


def mkdir(path):
    """Makes a directory if it does not already exist. Equivalent to bash `mkdir -p`."""
    if not os.path.exists(path):
        os.makedirs(path)
