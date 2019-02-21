#!/usr/bin/env python


# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import sys, os
import glob
import time
import datetime
import numpy as np
import astropy.io.fits as pyfits

import argparse
import pkg_resources
import appdirs

# Directory of this file
path = os.path.dirname(os.path.realpath(__file__))

# Global variables
import ac_settings as ac_set

# Location of ACTIN files:
actin_files_dir = os.path.join(path, "actin_files")
sys.path.append(actin_files_dir)
import ac_config
import ac_read_data
import ac_get_win
import ac_calc_ind
import ac_save
import ac_plot_time as ac_plot
import ac_tools

from matplotlib import pylab as plt #####

# initiate global variables:
ac_set.init()

# Configuration file:
config_file = os.path.join(path, "config_lines.txt")

# Version file:
version_file = os.path.join(path, "VERSION")

# Print preamble:
version = ac_set.preamble(version_file)


def actin_file(file, calc_index, rv_in, config_file=config_file, save_output=False, ln_plts=False, obj_name=None, targ_list=None, del_out=False, frac=True):
    """
    Runs ACTIN for one fits file.
    Accepts files of types: 'e2ds', 's1d', 's1d_*_rv', 'ADP', and 'rdb'.
    Recognizes fits files from HARPS and HARPS-N instruments.

    Parameters:
    -----------
    The same as the actin function below but for one file.

    Returns:
    --------
    data : dict
        Dictionary with data returned from fits files.

        Each key is a list with data related to a given measurement date
        given by the key 'date'.

        The used keys are:

        ==========  ========================================================
        keys        Description
        ----------  --------------------------------------------------------
        obj         str : Object (target) identification.
        flux        list : Flux per spectral order per pixel (e2ds) or flux
                    per pixel (other file types).
        wave        list : Wavelength per spectral order per pixel (e2ds) or
                    wavelength per pixel (other file types) [angstroms].
        blaze       list : Blaze function per spectral order (e2ds).
                    If blaze file is not found is list of 'ones'. For other
                    file types this header is not used.
        snr         {list, None} : SNR per spectral order. None if 'rdb'.
        median_snr  float : Median SNR of spectrum.
        date         str : Date of observation in the fits file format.
        bjd         float : Barycentric Julian date of observation [days].
        rv            float : Radial velocity [m/s] (if CCF file available).
        rv_err        float : Error on radial velocity (photon noise) [m/s]
                    (if CCF file available).
        fwhm        float : Full-Width-at-Half-Maximum of the CCF line
                    profile [m/s] (if CCF file available).
        cont        float : Contrast of the CCF line profile [%] (if CCF
                    file available).
        bis            float : Bisector Inverse Span of the CCF line profile
                    [m/s] (if BIS file available).
        noise        float : CCF noise [m/s] (if CCF file available).
        instr       str : Instrument identification.
        data_flg     str : Flag with value 'noDeblazed' when the blaze file
                    was not found (and flux_deb is real flux), None
                    otherwise.
        file_type   str : Type of file used as input: 'e2ds', 's1d'
                    (includes 's1d_*_rv'), 'ADP' or 'rdb'.
        ==========  ========================================================

    index : dict
        Dictionary containing the parameters related to the calculated spectral
        indices.

        Each key entry is a list of parameters where the list indices form the
        rows related to a given spectral index identified by the key 'index'.

        The returned keys are:

        ==========  ========================================================
        keys        Description
        ----------  --------------------------------------------------------
        index         str : Identification of the spectral index as given in
                    the configuration file.
        value         float : Spectral index value.
        error         float : Error of the index calculated by error propa-
                    gation.
        snr         {float, None} : Mean of the SNR at the lines spectral
                    order if the SNR per order was given as input, median
                    SNR of spectrum if median SNR given as input, 'None' if
                    no SNR values given as input.
        flg         {str, None} : Flags associated with the index: 'negFlux'
                    if negative flux detected in any of the lines used to
                    compute the index, None otherwise.
        mfrac_neg    list : Maximum fraction of flux with negative values
                    when taking into account all lines for a given index.
        ==========  ========================================================

    save_name : str
        Output rdb filename with path.
    """
    print()
    print("--------------------")
    print("EXECUTING ACTIN_FILE")
    print("--------------------")

    if type(file) is list: file = file[0]

    # Check if file is from known instrument
    tel, instr = ac_tools.get_instr(file)
    if instr == False: pass
    elif instr in ac_set.instr: pass
    else:
        msg="*** ERROR:\nUnrecognized instrument. ACTIN only accepts HARPS, HARPS-N or ESPRESSO. To read from a different instrument convert data to rdb file with the headers: 'obj', 'obs_date', 'bjd', 'wave', 'flux', 'error_pixel' (optional)."
        sys.exit(msg)

    # Checking if object in targ_list is the same as the object in fits file
    if targ_list:
        check = ac_tools.check_targ(file, targets=targ_list)
        if check is True: pass
        elif check is False: return

    # Read config file and retrieve lines information
    if calc_index:
        sel_lines = ac_config.read_conf(config_file, calc_index)

    # Read data from file
    data = ac_read_data.read_data(file, rv_in=rv_in, obj_name=obj_name)

    if not data: return

    # Check output file for duplicates
    if save_output is not False and data['file_type'] != "rdb":
        dupl = ac_save.check_duplicate(data['obj'], data['obs_date'], data['instr'], data['file_type'], save_output)
        if dupl is True: return

    info = {}
    info['config_file'] = config_file
    info['file_type'] = data['file_type']
    info['version'] = version
    info['source_path'] = os.path.split(file)[0]
    info['tel'] = data['tel']
    info['instr'] = data['instr']
    info['obj'] = data['obj']

    options = {}
    options['frac'] = frac

    if calc_index:
        # Check selected lines for spectral range and orders
        ac_calc_ind.check_lines(data['wave'], sel_lines)

        # Calculate flux in the required lines
        sel_lines = ac_calc_ind.calc_flux_lines(data, sel_lines, ln_plts=ln_plts, frac=frac)

        # Calculate chosen indices
        index = ac_calc_ind.calc_ind(sel_lines)
    if not calc_index:
        index = None
        sel_lines = None

    # Write output to rdb file in "out_dir"/"obj"
    if save_output is not False:
        rdb_file = ac_save.save_data(data, index, out_dir=save_output)
    else: rdb_file = None

    output = {}
    output['data'] = data
    output['index'] = index
    output['sel_lines'] = sel_lines
    output['info'] = info
    output['options'] = options
    output['rdb_file'] = rdb_file

    return output


def actin(files, calc_index, rv_in=None, config_file=None, save_output=False, ln_plts=False, obj_name=None, targ_list=None, del_out=False, frac=True, test=False):
    """
    Runs 'actin_file' function for one or multiple fits files, for one or multiple stars.
    Accepts fits files from HARPS, HARPS-N, and ESPRESSO instruments.
    Accepts files of types: 'S1D', 'S2D', 'e2ds', 's1d', 's1d_*_rv', 'ADP', and 'rdb'.
    Recognizes fits files from HARPS and HARPS-N instruments.
    """

    print()
    print("----------------")
    print(" STARTING ACTIN ")
    print("----------------")

    start_time = time.time()

    # Get config file from installation or input
    if config_file is None:
        cfg_file = get_config()
    else:
        cfg_file = config_file
    print()
    print("Using spectral lines from configuration file:")
    print(cfg_file)

    # test values can be 'S1D', 'S2D', 'e2ds', 's1d', 'adp', or 'rdb'
    # this needs to have before anything that uses 'files'
    if test:
        calc_index, files = ac_tools.test_actin(test, path)

    if not files:
        sys.exit()

    # Make lists to be iterated below
    files = list(files)

    if rv_in is None:
        rv_in = [rv_in]*len(files)
    else:
        rv_in = list(rv_in)

    # Check if files exist
    ac_tools.check_files(files)

    # Remove output file
    if del_out:
        ac_tools.remove_output(files, save_output, targ_list)

    # Option to make line plots directory the same as the data output dir
    if ln_plts == 'same':
        ln_plts = save_output

    total_files = len(files)

    # Organize files by path to star and file type
    lists_files = ac_tools.files_to_list_of_lists(files)

    for k in range(len(lists_files)):
        for i in range(len(lists_files[k])):
            n_files = 0
            for j in range(len(lists_files[k][i])):
                n_files += 1

                # Run actin file
                output = actin_file(lists_files[k][i][j],
                            calc_index,
                            rv_in=rv_in[j],
                            config_file=cfg_file,
                            save_output=save_output,
                            ln_plts=ln_plts,
                            obj_name=obj_name,
                            targ_list=targ_list,
                            del_out=del_out,
                            frac=frac)

            # POST-PROCESSING:
            if output:
                # Dictionaries of last file for file_type for each path
                sel_lines = output['sel_lines']
                info = output['info']
                options = output['options']
                rdb_file = output['rdb_file']

                # Save log and line info files
                ac_save.save_log(info, options, n_files, out_dir=save_output)
                ac_save.save_line_info(info, sel_lines, out_dir=save_output)

                # Save time-series plots
                ac_plot.plt_time(info, out_dir=save_output, rmv_flgs=False, save_plt=True)
                ac_plot.plt_time_mlty(info, out_dir=save_output, rmv_flgs=False, save_plt=True, hdrs=calc_index)
            else: pass

    elapsed_time = (time.time() - start_time)/60

    # Summary:
    print("\n---------------------------------")
    print("Fractional pixels:\t{}".format(frac))
    print("Files analysed:\t\t{}".format(total_files))
    print("Save output:\t\t{}".format(save_output))
    print("Elapsed time:\t\t{:.4f} min".format(elapsed_time))

    return


def get_config():
    """
    Check for existence of ACTIN folder and config file and creates them if not present. Returns the path to the config file.
    """
    cfg_dir = appdirs.user_config_dir('ACTIN')
    if not os.path.exists(cfg_dir):
        os.makedirs(cfg_dir)
    cfg_file = os.path.join(cfg_dir, 'config_lines.txt')
    if not os.path.isfile(cfg_file):
        create_user_config(cfg_file)

    return cfg_file


def create_user_config(cfg_file):
    """
    Create the user's config file
    """
    from shutil import copyfile ###
    src = pkg_resources.resource_stream(__name__, 'config_lines.txt')
    copyfile(src.name, cfg_file)



def main():
    """
    Main function, call actin function with arguments from terminal.
    """

    # initiate the parser
    parser = argparse.ArgumentParser()

    # add short and long argument
    parser.add_argument('--files', '-f', help='Read file(s)', nargs='+')

    parser.add_argument('--calc_index', '-i', help="Index id to calculate as designated by 'ind_id' in config_index.txt.", nargs='+', default=None)

    parser.add_argument('--rv_in', '-rv', help="RV value to calibrate wavelength. If False (default) try to read RV from CCF file.", nargs='+', default=None, type=float)

    parser.add_argument('--config_file', '-cf', help='Path to config_file, or False (default) read config file from standard directory.', default=None)

    parser.add_argument('--save_output', '-s', help='Path to output directory of data table, or False (default).', default=False)

    parser.add_argument('--ln_plts', '-lp', help="Path to directory to save line plots. If 'same' saves line plots to same directory of data output. If 'show' only shows the plots. If 'False' (default) does not save or show line plots", default=False, type=lambda x: (str(x).lower() == 'true'))

    parser.add_argument('--obj_name', '-obj', help='Give target a name that overides the one from the fits files.', default=None)

    parser.add_argument('--targ_list', '-tl', help='Give a list of stars to select from fits files.', nargs='+', default=None)

    parser.add_argument('--del_out', '-del', help='Delete output data file if True.', default=False, type=lambda x: (str(x).lower() == 'true'))

    parser.add_argument('--test', '-t', help='Tests actin using the provided fits files in the "test_files" directory. Options are "e2ds", "s1d", and "adp"', default=False, type=lambda x: (str(x).lower() == 'true'))

    parser.add_argument('--frac', '-frc', help='Turns fractional pixel on (True, default) or off (False).', default=True, type=lambda x: (str(x).lower() == 'true'))

    #parser.add_argument('--plt_spec', '-pspec', help='Plot full spectrum if True, If int is given plot spectrum in the int order. False is default.', default=False)

    # read arguments from the command lines
    args = parser.parse_args()

    actin(files=args.files,
        calc_index=args.calc_index,
        rv_in=args.rv_in,
        config_file=args.config_file,
        save_output=args.save_output,
        ln_plts=args.ln_plts,
        obj_name=args.obj_name,
        targ_list=args.targ_list,
        del_out=args.del_out,
        test=args.test,
        frac=args.frac)

if __name__ == "__main__":
    main()
