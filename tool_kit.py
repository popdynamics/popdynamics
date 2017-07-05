
import platform
import os
import glob

"""
General static methods for use in any example module.
"""


def ensure_out_dir(out_dir):
    """
    Make sure the output directory exists and create if it doesn't.
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)


def open_out_dir(out_dir):
    """
    Open the output dir at the end of the model run.
    """

    pngs = glob.glob(os.path.join(out_dir, '*png'))
    operating_system = platform.system()
    if 'Windows' in operating_system:
        os.system("start " + " ".join(pngs))
    elif 'Darwin' in operating_system:
        os.system('open ' + " ".join(pngs))