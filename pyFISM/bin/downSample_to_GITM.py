#/usr/bin/env python

import os
import re
import sys
import numpy as np

import argparse
import datetime
import glob

import scipy.io

import pycurl
import configparser


_rpath_ = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.insert(1, os.path.join(_rpath_, 'lib'))


import utils as utl



def get_FISM_results(archive, ydoy, remURL, verb):

    syear = str(ydoy.year)
    strydoy = ydoy.strftime("%Y%j")

    savfiles = glob.glob(os.path.join(archive, syear, 'FISM_60sec_'+strydoy+'*.sav'))
    savfiles.sort()
    nsav = len(savfiles)
    savefile = None
    if savfiles is None or nsav == 0:
        if verb:
            print(f"Retrieving file from remote site {remURL}")
        savefile = utl.retrieveSAVfile(archive, ydoy, remURL, verb)
    else:
        savefile = savfiles[-1]
        if verb:
            print(f"Using existing file {savefile}")

    if not os.path.exists(savefile):
        print(f"ERROR - Unable to find or download file {savefile}")
        sys.exit(1)

    savdata = scipy.io.readsav(savefile, python_dict=True)

    savdata['savefile'] = savefile

    return savdata



def ConvertPar(ConDict, parName, type='float', verb=0):
    if verb > 1:
        print(f"Converting parameter {parName} to a list of floats")
    lines = ConDict[parName].split("\n")
    rList = []
    for line in lines:
        items = re.split(r'\s+', line)
        for val in items:
            if type == 'float':
                rList.append(float(val))
            elif type == 'int':
                rList.append(int(val))
            elif type == 'str':
                rList.append(val)

    return rList





#FISM_60sec_2023001_v02_01.sav

def main(args):
    yd = args.yd
    verb = args .verbose
    archive = args.archive
    config = args.config
    bin_file = args.bin_file



    if not os.path.exists(config):
        print(f"ERROR - Unable to find default config file {config}")
        sys.exit(1)

    Config = configparser.ConfigParser()
    Config.optionxform = str
    Config.read(config)

    ConDict = (utl.ConvertDict(Config, verb))['defaults']
    if archive == '':
        archive = ConDict['archive']


    rootdir = ConDict['rootdir'] #'/Users/romeog1/BOX/eclipse/'
    rootURL = ConDict['rootURL']



    # Restore the found best primary proxy vs wavelength
    # EUV/EVE
    best_primary_proxy = ConDict['best_primary_proxy']
    best_primary_proxy_fuv = ConDict['best_primary_proxy_fuv']
    bestpp = scipy.io.readsav(os.path.join(archive, best_primary_proxy), python_dict=True)
    betsFUV = scipy.io.readsav(os.path.join(archive, best_primary_proxy_fuv), python_dict=True)


    start_wv = ConvertPar(ConDict, 'start_wv', 'float', verb)
    end_wv = ConvertPar(ConDict, 'end_wv', 'float', verb)

    ngitmbins = len(start_wv)

    names = ConvertPar(ConDict, 'paramNames', 'str', verb)
    n_names = len(names)

    if not os.path.exists(archive):
        try:
            os.makedirs(archive, exist_ok=True)
        except:
            print(f"ERROR - Unable to create directory archive {archive}")
            sys.exit(1)

    # Create a best aia mask tag for all FISM wavelenghts
    # 0:94; 1:171, 2:193, 3:211, 4:304, 5:335, 6:1600
    msk_tg = np.zeros(1900, dtype=np.int32)
    #   ; Set XUV to the 94A (tag=0), or XRT (tag=7) if/when avaialble
    msk_tg[0:59] = 7
    # Adjust to match EVE best proxy
    # 1'mgii' -> 304
    # 2'f107' -> 94
    # 3'goes' -> XRT
    # 4'lya' -> 304
    # 5'QD' -> XRT
    # 6'171' -> 171
    # 7'304' -> 304
    # 8'335' -> 335
    # 9'369' -> 335
    # 10'171d' -> 171
    # 11'304d' -> 304
    # 12'lyad' -> 304
    for i in range(60, 1060):
        match bestpp['best_primary_tag'][i - 60]:
            case 1:
                msk_tg[i] = 4
            case 2:
                msk_tg[i] = 0
            case 3:
                msk_tg[i] = 7
            case 4:
                msk_tg[i] = 4
            case 5:
                msk_tg[i] = 7
            case 6:
                msk_tg[i] = 1
            case 7:
                msk_tg[i] = 4
            case 8:
                msk_tg[i] = 5
            case 9:
                msk_tg[i] = 5
            case 10:
                msk_tg[i] = 1
            case 11:
                msk_tg[i] = 4
            case 12:
                msk_tg[i] = 4

    # Set all FUV to 304A, 1600A, or 1700A for now
    msk_tg[1060:1400] = 4
    msk_tg[1401:1650] = 6
    msk_tg[1651:1899] = 9

    # Force 13-14nm to be 131A and 28-29nm to be 284A
    msk_tg[130:140] = 8
    msk_tg[280:290] = 10




    nyd = len(yd)
    ydoy = None
    if nyd == 7:
        ydoy = datetime.datetime.strptime(yd, "%Y%j")
    elif nyd == 8:
        ydoy = datetime.datetime.strptime(yd, "%Y%m%d")
    else:
        print(f"ERROR - Unable to parse input date {yd}")
        sys.exit(1)



    strymd = ydoy.strftime("%Y%m%d")
    syear = str(ydoy.year)
    strydoy = ydoy.strftime("%Y%j")
    sdoy = ydoy.strftime("%j")

    # Read (or first get) IDL saveset file
    fism_res = get_FISM_results(archive, ydoy, rootURL, verb)


    ##########################################################################
    #
    # Doing TLSM downsampling


    # read input solar spectra.
    spectra1 = utl.get_spectra_saveset(fism_res, verb)

    if spectra1 is None:
        sys.exit(1)

    # scale to 1 nm binning scheme
    spectra2 = utl.hi_to_nm(spectra1, _rpath_, verb)
    if spectra2 is None:
        sys.exit(1)

    # scale to Hinteregger resolution
    in_spectra = utl.nm_to_hs(spectra2, _rpath_, verb)

    # rebin based on selected binning scheme.
    out_spectra = utl.rebin(bin_file, in_spectra, _rpath_, verb)

    # output
    utl.write_spectra_saveset(out_spectra, fism_res, 'TLSM', verb)

    ##########################################################################
    #
    # Doing GITM downsampling

    irradiance = fism_res['irradiance']
    irr_shape = list(irradiance.shape)
    nmin = irr_shape[0]

    msk_irradiance_gitm = np.zeros([nmin, 59], dtype=np.float32)  # declare array for GITM binned data

    w_gitm_list = []
    for m in range(ngitmbins):
        w_gitm_list.append(np.where((fism_res['wavelength'] >= start_wv[m]) & (fism_res['wavelength'] < end_wv[m]))[0])

    for j in range(0, nmin):
        # Put full masked spectrum into GITM bins
        for m in range(ngitmbins):
            ir_gitm = np.sum(irradiance[j,w_gitm_list[m]])/10.    #/10 as FISM2 has 0.1nm bins
            msk_irradiance_gitm[j,m] = ir_gitm            # W/m^2


    gitm_out_spectra = np.zeros([msk_irradiance_gitm.shape[0]+2, msk_irradiance_gitm.shape[1]], dtype=np.float64)
    gitm_out_spectra[0,:] = start_wv
    gitm_out_spectra[1,:] = end_wv
    gitm_out_spectra[2:,:] = msk_irradiance_gitm[:,:]

    utl.write_spectra_saveset(gitm_out_spectra, fism_res, 'GITM', verb)

    pass




if __name__ == "__main__":

    #fism_for_eclipse, yd = yd, msk_fl = msk_fl, fl_time = fl_time
    parser = argparse.ArgumentParser(
        prog="fism_for_eclipse", description="Python version of the IDL code"
    )

    parser.add_argument(
        "-D","--yd",
        help="""YearDOY in the for YYYYDDD (e.g. 2023287)""",
        #default='2023287',
        default='2024097',
        type = str
    )


    parser.add_argument(
        "-f","--bin_file",
        help="""Select binning scheme""",
        default=os.path.join(_rpath_,'input','low_bins.txt'),
        type = str
    )



    parser.add_argument(
        "-C", "--config",
        help="""Configuration file""",
        required=False,
        default=os.path.join(_rpath_, 'config', 'fism_for_eclipse.config'),
        type = str
    )

    parser.add_argument(
        "-A", "--archive",
        help="""Archive of downloaded mask files""",
        required=False,
        #default='/Users/romeog1/BOX/eclipse/FISM_Archive',
        default='',
        type = str
    )

    parser.add_argument(
        "-V", "--verbose",
        help="""Verbosity level, default is 0""",
        required=False,
        default=0,
        type = int
    )


    args = parser.parse_args()

    main(args)