#!/usr/bin/env python

# compatibility with python 2/3:
from __future__ import print_function
from __future__ import division

import os, sys
import string
from pylab import *
import astropy.io.fits as pyfits
import numpy as np

from matplotlib import pylab as plt

import ac_settings as ac_set


def compute_flux(wave, flux, blaze, noise, ln_ctr, ln_win, bandtype, frac=True, test_plot=False):

    print("Executing: compute_flux")

    if blaze is None: blaze = np.ones(len(flux))

    # make all important values the size of px_size
    px_size = np.diff(wave)
    #px_size = [float(px_size[k]) for k in range(len(px_size))] ####
    #px_size = np.asarray(px_size) ####
    wave = wave[1:]
    flux = flux[1:]
    blaze = blaze[1:]

    # make all values float otherwise python will not detect very small differences between wave and wmin/max
    #wave = [float(wave[k]) for k in range(len(wave))] ###
    #wave = np.asarray(wave) ####

    # BANDPASS TYPES
    if bandtype == 'tri':
        #wmin = float(ln_ctr - ln_win)
        #wmax = float(ln_ctr + ln_win)
        wmin = ln_ctr - ln_win #########
        wmax = ln_ctr + ln_win #########
        cond = (wave > wmin) & (wave < wmax)
        bandfunc = -np.abs(wave-ln_ctr)/ln_win + 1.
        bandfunc = np.where(bandfunc > 0, bandfunc, bandfunc*0.0)
    if bandtype == 'sq':
        #wmin = float(ln_ctr - ln_win/2.)
        #wmax = float(ln_ctr + ln_win/2.)
        wmin = ln_ctr - ln_win/2. #########
        wmax = ln_ctr + ln_win/2. #########
        cond = (wave > wmin) & (wave < wmax)
        bandfunc = np.ones(len(wave))
        bandfunc = np.where(cond, 1., 0.)

    # HARPS METHOD
    if frac == False:
        flux_win = flux[cond]
        wave_win = wave[cond]
        blaze_win = blaze[cond]
        flux_win_deb = flux[cond]/blaze[cond]
        px_size_win = px_size[cond]
        npixels = len(wave[cond])
        response = bandfunc[cond]
    # ACTIN METHOD
    if frac == True:
        dwave_l = wave[cond][0]-wmin
        dwave_r = wmax-wave[cond][-1]
        dwave = dwave_l + dwave_r

        cond_l = (wave < wave[cond][0])
        cond_r = (wave > wave[cond][-1])

        frac_l = dwave_l/px_size[cond_l][-1]
        frac_r = dwave_r/px_size[cond_r][0]

        npixels = len(wave[cond]) + frac_l + frac_r

        # just for visualization purposes
        wave_win = np.r_[wave[cond][0]-dwave_l, wave[cond], wave[cond][-1]+dwave_r]

        flux_win = np.r_[flux[cond_l][-1], flux[cond], flux[cond_r][0]]

        blaze_win = np.r_[blaze[cond_l][-1], blaze[cond], blaze[cond_r][0]]

        px_size_win = np.r_[px_size[cond_l][-1], px_size[cond], px_size[cond_r][0]]

        response = bandfunc[cond]
        if bandtype == 'sq':
            response = np.r_[1*frac_l, response, 1*frac_r]
        if bandtype == 'tri':
            response = np.r_[0.0, response, 0.0]

        flux_win_deb = flux_win/blaze_win

        # tests
        bandwidth = dwave_l + wave[cond][-1]-wave[cond][0] + dwave_r
        if bandwidth == (wmax - wmin): pass
        else: print("*** bandwidth test 1 ERROR")

        bandwidth = frac_r*px_size_win[-1] + sum(px_size[cond][1:]) + px_size_win[0]*frac_l
        if bandwidth == (wmax - wmin): pass
        else: print("*** bandwidth test 2 ERROR")

    # Flux sum and variance for line:
    f_sum = sum(flux_win_deb*response/sum(response*px_size_win))
    f_sum_var = sum((flux_win+noise**2)*response**2/blaze_win**2)/sum(response*px_size_win)**2

    print("Pixels in bandpass = %.4f" % npixels)

    # Flag negative flux inside bandpass
    flg, frac_neg = flag_negflux(flux_win)

    # Test plots
    if test_plot == True:
        w_out_l = wave[wave < wave[cond][0]][-1]
        w_in_l = wave[cond][0]
        w_out_r = wave[wave > wave[cond][-1]][0]
        w_in_r = wave[cond][-1]

        plt.plot(wave, flux, 'kx')
        plt.plot(wave_win, flux_win*response, 'b.')
        plt.plot(wave_win, response*max(flux_win),c='orange',ls='-',linewidth=2)

        plt.axhline(0.0, c='k', ls=':')

        plt.axvline(wmin, color='k', ls='-')
        plt.axvline(w_in_l, color='r', ls=':')
        plt.axvline(w_out_l, color='g', ls='--')

        plt.axvline(wmax, color='k', ls='-')
        plt.axvline(w_in_r, color='r', ls=':')
        plt.axvline(w_out_r, color='g', ls='--')
        plt.show()

    return f_sum, f_sum_var, bandfunc, npixels, flg, frac_neg


def flag_negflux(flux):
    """
    Tests if flux has negative values and returns flag 'flg' as 'negFlux' if detected, None otherwise, and the fraction of pixels with negative values of flux, 'frac_neg'.
    """
    negflux_array = np.where(flux < 0.0, flux, 0.0)
    #negflux = sum(negflux_array)

    negflux_only = [negflux_array[x] for x in range(len(negflux_array)) if negflux_array[x] < 0.0]

    # fraction of pixels with negative flux
    frac_neg = len(negflux_only)/len(flux)

    #frac_neg = sum(abs(negflux))/sum(flux)

    flag_array = np.where(flux < 0.0, 'negFlux', None)

    if 'negFlux' in flag_array:
        flg = 'negFlux'
        print("*** WARNING: Negative flux detected")
        print("Fraction of pixels with negative values = {:.5f}".format(frac_neg))
    else: flg = None

    return flg, frac_neg


def remove_output(files, save_output, targ_list=None):
    """
    Removes output directories for given fits files 'files', 'save_output' directory and list of targets 'targ_list' (if available).
    """

    print()
    print("Executing actin_functions.remove_output:")
    print("Searching output files to delete...")

    fn_rdb = ac_set.fnames['data']
    fn_log = ac_set.fnames['log_data']
    fn_lns = ac_set.fnames['lines_data']
    ftypes = ac_set.ftypes['all']

    instr_list = []
    for k in range(len(files)):
        instr_list.append(get_instr(files[k])[1])
    instr = list(set(instr_list))

    if isinstance(files,str): files = [files]

    err_msg = None

    if save_output is False:
        print("*** ERROR: No 'save_output' directory to search files.\nInsert the path to output folder after '-s'.")
        sys.exit()

    obj = []
    ft = []
    for file in files:
        for type in ftypes:
            if type in file:
                ft.append(type)
                if type == "rdb":
                    obj.append(read_rdb(file)[0]['obj'][0])
                elif type != "rdb":
                    obj.append(get_target(file))

    if targ_list is not None:
        ft = [ft[k] for k in range(len(ft)) if obj[k] in targ_list]
        obj = [obj[k] for k in range(len(ft)) if obj[k] in targ_list]

    # remove rdb and log files
    old_obj = None
    old_ft = None
    for k in range(len(obj)):
        for i in range(len(instr)):
            if obj[k] != old_obj or ft[k] != old_ft:
                print("Object: {}".format(obj[k]))
                print("File type: {}".format(ft[k]))
                file_rdb = "{}_{}_{}_{}".format(obj[k], instr[i], ft[k], fn_rdb)
                file_log = "{}_{}_{}_{}".format(obj[k], instr[i], ft[k], fn_log)
                file_lns = "{}_{}_{}_{}".format(obj[k], instr[i], ft[k], fn_lns)
                ff_rdb = os.path.join(save_output, obj[k], file_rdb)
                ff_log = os.path.join(save_output, obj[k], file_log)
                ff_lns = os.path.join(save_output, obj[k], file_lns)

                if os.path.isfile(ff_rdb):
                    os.remove(ff_rdb)
                    print("Output file removed: {}".format(ff_rdb))
                if os.path.isfile(ff_log):
                    os.remove(ff_log)
                    print("Output file removed: {}".format(ff_log))
                if os.path.isfile(ff_lns):
                    os.remove(ff_lns)
                    print("Output file removed: {}".format(ff_lns))

                else:
                    print("No files to remove.")
            else: pass
            old_obj = obj[k]
            old_ft = ft[k]

    return


def files_to_list_of_lists(files):
    """
    Organize files by path and file type.
    """
    # path for each star in files
    source_path = []
    for k in range(len(files)):
        source_path.append(os.path.split(files[k])[0])
    source_path = list(set(source_path))

    # different file types in files
    ft = []
    for k in range(len(files)):
        ft.append(get_file_type(files[k]))
    ft = list(set(ft))

    # organize files by path and file type
    lists_files = []
    for k in range(len(source_path)):
        ft_files = []
        for k in range(len(ft)):
            ft_f = []
            for i in range(len(files)):
                ftype = get_file_type(files[i])
                if ftype == ft[k]:
                    ft_f.append(files[i])
            ft_files.append(ft_f)
        lists_files.append(ft_files)

    return lists_files


def override_obj(obj, obj_name):
    # Override object name with name given in obj_name option
    if type(obj_name) is list and len(obj_name) == 1:
        return obj_name[0]
    elif type(obj_name) is list and len(obj_name) > 1:
        print("*** ERROR: obj_name requires only one name, more than one given")
        return
    else: return obj_name


def check_files(files):
    for k in range(len(files)):
        if not os.path.isfile(files[k]):
            print("*** ERROR: File {} not found.".format(files[k]))
            sys.exit()


def get_file_type(file):
    """
    Checks file name for known file types and returns it.
    """
    for k in range(len(ac_set.ftypes['all'])):
        if ac_set.ftypes['all'][k] in file:
            return ac_set.ftypes['all'][k]


def get_instr(fits_file):
    if ".rdb" not in fits_file:
        hdu = pyfits.open(fits_file)

        tel = hdu[0].header['TELESCOP']
        instr = hdu[0].header['INSTRUME']

        hdu.close()

        return tel, instr
    if ".rdb" in fits_file: return False, False


def get_target(fits_file):
    """
    Returns the object targeted in the fits file 'fits_file'.
    """

    hdu = pyfits.open(fits_file)

    try:
        obj = hdu[0].header['OBJECT']
        hdu.close()
    except:
        try:
            obj = hdu[0].header['ESO OBS TARG NAME']
            hdu.close()
        except:
            try:
                obj = hdu[0].header['TNG OBS TARG NAME']
                hdu.close()
            except:
                print("*** ERROR: Cannot identify object.")
                return

    return obj


def check_targ(fits_file, targets):
    """
    Checks if a fits file belongs to a target in a list of targets.
    """

    print("Executing: check_targ")

    print("Targets = {}".format(targets))

    obj = get_target(fits_file)

    print("Object = {}".format(obj))
    if obj in targets: return True
    else:
        print("*** INFO: {} is not in the target list.".format(obj))
        print("*** INFO: file {}".format(fits_file))
        return False


def test_actin(test, path):
    calc_index = ("I_CaII", "I_Ha", "I_NaI")
    if test == "S1D":
        files = os.path.join(path, "test_files", "HD41248_1_1_S1D_A.fits")
    elif test == "S2D":
        files = os.path.join(path, "test_files", "HD41248_1_1_S2D_A.fits")
    elif test == "e2ds":
        files = os.path.join(path, "test_files", "HARPS.2003-12-13T06:19:48.371_e2ds_A.fits")
    elif test == "s1d":
        files = os.path.join(path, "test_files", "HARPS.2010-09-18T23:42:36.178_s1d_A.fits")
    elif test == "adp":
        files = os.path.join(path, "test_files", "ADP.2014-09-16T11:04:45.123.fits")
    elif test == 'rdb':
        files = os.path.join(path, "test_files", "2010-09-18T23:42:36.178_spec.rdb")
    else:
        print("*** ERROR:")
        print("*** Test can only be 'S1D', 'S2D', 'e2ds', 's1d', 'adp', or 'rdb'")
    return calc_index, files


def read_rdb(filename):
    """
    Reads an .rdb file and returns a dictionary with the headers as keys and data as lists ('output'), and also a list of headers ('keys').
    """
    # use: table = pyrdb.read_rdb(file)[0] for data
    # use: table = pyrdb.read_rdb(file)[1] to get the keys by order

    f = open(filename, 'r')
    data = f.readlines()
    f.close()

    key = str.split(data[0][:-1],'\t')
    output = {}
    for i in range(len(key)): output[key[i]] = []

    for line in data[2:]:
        qq = str.split(line[:-1],'\t')
        for i in range(len(key)):
            try: value = float(qq[i])
            except ValueError: value = qq[i]
            output[key[i]].append(value)

    return output, key


def save_rdb(dic, keys,file):
    """
    From a disctionary 'dic' saves the columns related to the specified 'keys' into an .rdb file named 'file'.
    """

    out = open(file,'w')
    n_keys = len(keys)

    for k in range(n_keys):
        if k != n_keys-1: out.write('{}\t'.format(keys[k]))
        else: out.write('{}'.format(keys[k]))
    out.write('\n')
    for k in range(n_keys):
        if k != n_keys-1: out.write('{}\t'.format('-'*len(keys[k])))
        else: out.write('{}'.format('-'*len(keys[k])))
    out.write('\n')

    for i in range(len(dic[keys[0]])):
        for k in range(n_keys):
            if k != n_keys-1: out.write('{}\t'.format(dic[keys[k]][i]))
            else: out.write('{}'.format(dic[keys[k]][i]))
        out.write('\n')
    out.close()


def add_rdb(dic,keys, file_name):
    """
    Adds data to an existing .rdb file. The 'keys' must match the headers already present in the file.
    """

    out = open(file_name,'a')
    n_keys = len(keys)
    for i in range(len(dic[keys[0]])):
        for k in range(n_keys):
            if k != n_keys-1: out.write('{}\t'.format(dic[keys[k]][i]))
            else: out.write('{}\t'.format(dic[keys[k]][i]))
        out.write('\n')
    out.close()


def plot_params(width=6, height=3.5):
        """
        Parameters for plots.
        """

        rcdefaults()
        rc('axes', linewidth=1)
        rc('lines', markeredgewidth=0.5)

        rc('errorbar', capsize=2)

        rc('ytick', direction='in')
        rc('xtick', direction='in')

        rc('ytick', right='True')
        rc('xtick', top='True')

        rc('xtick.minor', visible='True')


        rc('xtick.major', size=7)
        rc('ytick.major', size=7)
        rc('xtick.minor', size=3)
        rc('ytick.minor', size=3)
        rc("font", family="sans-serif")
        rc("font", size=10)

        rc("figure.subplot", left=(0.15))
        rc("figure.subplot", right=(0.93))
        rc("figure.subplot", bottom=(0.12))
        rc("figure.subplot", top=(0.95))

        return width, height
