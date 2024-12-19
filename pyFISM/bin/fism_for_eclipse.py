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
    mask_fl = args.msk_fl
    fl_time = args.fl_time
    verb = args .verbose
    archive = args.archive
    config = args.config

    if config is None or config == '' or not os.path.exists(config):
        _rpath_ = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
        fname = (os.path.basename(__file__)).replace('.py', '.config')
        config = os.path.join(_rpath_, 'config', fname)
        if verb:
            print(f"INFO - Trying to use default confi file {config}")

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

    fl_times = []
    if fl_time == 'all':
        cydoy = ydoy
        nydoy = ydoy + datetime.timedelta(days=1)
        while cydoy <= nydoy:
            fl_times.append(cydoy.strftime("%H%M"))
            cydoy += datetime.timedelta(minutes=5)
    else:
        fl_times.append(fl_time)



    strymd = ydoy.strftime("%Y%m%d")
    syear = str(ydoy.year)
    strydoy = ydoy.strftime("%Y%j")
    sdoy = ydoy.strftime("%j")

    fism_res = get_FISM_results(archive, ydoy, rootURL, verb)
    for fltime in fl_times:
        if mask_fl == '' or mask_fl is None or fl_time == 'all':
            mask_fl = os.path.join(rootdir, 'Eclipse_'+strymd, 'mask', 'fism2-input', strymd+fltime+'00_150km_alleof.nc')

        if not os.path.exists(mask_fl):
            print(f"ERROR - Unable to find input file {mask_fl}")
            continue


        pDict = {}
        pDict['netfile'] = mask_fl
        pDict['ValuesOnly'] = 0
        pDict['getMetadata'] = 1
        pDict['reMapNames'] = 0
        jnk, nc_data = utl.ReadNetcdfFile(pDict)

        fname = os.path.basename(mask_fl)
        msk_st_hr = int(fname[8:10])
        msk_st_min = int(fname[10:12])
        msk_st_sec = int(fname[12:14])
        msk_st_utc = msk_st_hr * 3600 + msk_st_min * 60 + msk_st_sec

        # Find the numbers of lats and longs in the file
        nlat = nc_data['dimensions']['glat']
        nlon = nc_data['dimensions']['glon']


        fism_st_ind_tmp = np.where(fism_res['utc'] >= msk_st_utc)[0]
        fism_st_ind = fism_st_ind_tmp[0]
        sod = fism_res['utc'][fism_st_ind_tmp[0]]

        fism2 = fism_res['irradiance'][fism_st_ind, :]

        sza = np.zeros([nlat, nlon], dtype=np.float32)  # declare solar zenith angle array
        iseclipse = np.zeros(sza.shape)



        # Apply the mask factors vs wavelength to the FISM spectra
        msk_irradiance = np.zeros([nlat, nlon, 1900], dtype=np.float32)     # declare array
        msk_irradiance_gitm = np.zeros([nlat, nlon, 59], dtype=np.float32)  # declare array for GITM binned data

        w_gitm_list = []
        for m in range(ngitmbins):
            w_gitm_list.append(np.where((fism_res['wavelength'] >= start_wv[m]) & (fism_res['wavelength'] < end_wv[m]))[0])

        for j in range(0, nlat):
            for k in range(0, nlon):
                msk_fct = np.zeros(1900, dtype=np.float32)
                #names = ['94', '171', '193', '211', '304', '335', '1600', 'X', '131', '1700', '284']
                # for i in range(n_names):
                #     idxs = np.where(msk_tg == i)[0]
                #     nidx = len(idxs)
                #     if nidx > 0:
                #         msk_fct[idxs] = nc_data['variables'][names[i]]['values'][j, k]

                if sza[j,k] < 90.0:
                    a = np.ones(n_names)
                    for i in range(n_names):
                        if names[i] not in nc_data['variables']:
                            continue

                        a[i] = nc_data['variables'][names[i]]['values'][j,k]

                        idxs = np.where(msk_tg == i)[0]
                        nidx = len(idxs)
                        if nidx > 0:
                            msk_fct[idxs] = nc_data['variables'][names[i]]['values'][j, k]

                    if not np.logical_or(np.all(a==0), np.all(a==1)):
                        iseclipse[j,k] = 1

                msk_irradiance[j, k, :] = fism_res['irradiance'][fism_st_ind, :] * msk_fct
                sza[j,k] = nc_data['variables']['sza']['values'][j,k]

                # Put full masked spectrum into GITM bins
                for m in range(ngitmbins):
                    #w_gitm = np.where((fism_res['wavelength'] >= start_wv[m]) & (fism_res['wavelength'] < end_wv[m]))[0]
                    ir_gitm = np.sum(msk_irradiance[j,k,w_gitm_list[m]])/10.    #/10 as FISM2 has 0.1nm bins
                    fism2_gitm = np.sum(fism2[w_gitm_list[m]])/10    #/10 as FISM2 has 0.1nm bins Convert input fism2 spectra
                    msk_irradiance_gitm[j,k,m] = ir_gitm            # W/m^2

                # for i in range(0, 1899):
                #     match msk_tg[i]:
                #         case 0: msk_fct = nc_data['variables']['94']['values'][j, k]
                #         case 1: msk_fct = nc_data['variables']['171']['values'][j, k]
                #         case 2: msk_fct = nc_data['variables']['193']['values'][j, k]
                #         case 3: msk_fct = nc_data['variables']['211']['values'][j, k]
                #         case 4: msk_fct = nc_data['variables']['304']['values'][j, k]
                #         case 5: msk_fct = nc_data['variables']['335']['values'][j, k]
                #         case 6: msk_fct = nc_data['variables']['1600']['values'][j, k]
                #         case 7: msk_fct = nc_data['variables']['X']['values'][j, k]
                #         case 8: msk_fct = nc_data['variables']['131']['values'][j, k]
                #         case 9: msk_fct = nc_data['variables']['1700']['values'][j, k]
                #         case 10: msk_fct = nc_data['variables']['284']['values'][j, k]

                #    msk_irradiance[j, k, i] = fism_res['irradiance'][fism_st_ind, i] * msk_fct



        # Save data into arrays
        lat = nc_data['variables']['glat']['values']
        lon = nc_data['variables']['glon']['values']

        outdir = os.path.join(archive, 'aia_masks', syear, sdoy+'_new2')
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)

        outfile1 = os.path.join(outdir, 'FISM_MASKED_'+strymd+'_'+fltime+'00.nc4')
        outfile2 = os.path.join(outdir, 'FISM_MASKED_'+strymd+'_'+fltime+'00_GITM.nc4')

        # dict1 = {'msk_irradiance':msk_irradiance, 'sod':sod, 'ydoy':strydoy,
        #          'lat':lat, 'lon':lon, 'wavelength':fism_res['wavelength']}
        #
        # dict2 = {'msk_irradiance_gitm':msk_irradiance_gitm, 'sod':sod, 'ydoy':strydoy,
        #          'lat':lat, 'lon':lon, 'start_wv':start_wv, 'end_wv':end_wv}


        dict1 = {'msk_irradiance':msk_irradiance, 'sod':sod, 'ydoy':strydoy,
                 'lat':lat, 'lon':lon, 'wavelength':fism_res['wavelength'],
                 'sza':sza, 'iseclipse':iseclipse, 'fism2':fism2}

        dict2 = {'msk_irradiance_gitm':msk_irradiance_gitm, 'sod':sod, 'ydoy':strydoy,
                 'lat':lat, 'lon':lon, 'start_wv':start_wv, 'end_wv':end_wv,
                 'sza':sza, 'iseclipse':iseclipse, 'fism2': fism2, 'wavelength':fism_res['wavelength']}



        compDict = None
        if 'compare' in ConDict and ConDict['compare'] != '0':
            # compfile = os.path.join(rootdir, 'Eclipse_20240408', 'mask', 'fism2-output', 'fism_masked_'+\
            #                         strymd+'_'+fltime+'00.sav')
            compfile = os.path.join(archive, 'aia_masks', 'idl', syear, 'fism_masked_'+\
                                    strymd+'_'+fltime+'00.sav')
            if os.path.exists(compfile):
                compDict = scipy.io.readsav(compfile, python_dict=True)
                newar = np.transpose(compDict['msk_irradiance'])

                try:
                    diff = (newar - dict1['msk_irradiance'])/newar
                except:
                    print(f"ERROR - Unable to calculate diff")



        utl.writeFISMNetCDF(outfile1, dict1, verb)
        utl.writeFISMNetCDF(outfile2, dict2, verb)

        pass


if __name__ == "__main__":

    #fism_for_eclipse, yd = yd, msk_fl = msk_fl, fl_time = fl_time
    parser = argparse.ArgumentParser(
        prog="fism_for_eclipse", description="Python version of the IDL code"
    )

    parser.add_argument(
        "--yd",
        help="""YearDOY in the for YYYYDDD (e.g. 2023287)""",
        #default='2023287',
        default='2024099',
        type = str
    )

    parser.add_argument(
        "--msk_fl",
        help="""Input netCDF file.""",
        default = '',
        type = str
    )

    parser.add_argument(
        "-T", "--fl_time",
        help="""Event time in the form of HHMM (e.g. 0650), Using 'all' will make the script do all times, with increments of 5 minutes""",
        required=False,
        default='1920',
        type = str
    )

    parser.add_argument(
        "-C", "--config",
        help="""Configuration file""",
        required=False,
        default='',
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

# -V 1 -T 1920

