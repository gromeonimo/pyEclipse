#!/usr/bin/env python3

import os
import sys

import argparse
import datetime


_rpath_ = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.insert(1, os.path.join(_rpath_, 'lib'))


import utils as utl
import glob





def main(args):
    bin_file = args.bin_file
    indir = args.indir
    filepatt = args.filepatt
    verb = args.verbose

    infiles = glob.glob(os.path.join(indir, filepatt))
    infiles.sort()

    for solar_file in infiles:
        if verb:
            print(f"INFO - Processing file {solar_file}")
    # read input solar spectra.
        spectra1, nc_data = utl.get_spectra_fism2(solar_file, verb)

        if spectra1 is None:
            sys.exit(1)

    # scale to 1 nm binning scheme
        spectra2 = utl.hi_to_nm(spectra1, _rpath_, verb)
        if spectra2 is None:
            sys.exit(1)

    # scale to Hinteregger resolution
        in_spectra = utl.nm_to_hs(spectra2, _rpath_, verb)

    # rebin based on selected binning scheme.
        out_spectra = utl.rebin(bin_file,in_spectra, _rpath_, verb)



    # output
        utl.put_spectra_fism2(solar_file, out_spectra, nc_data, verb)



if __name__ == "__main__":


    #fism_for_eclipse, yd = yd, msk_fl = msk_fl, fl_time = fl_time
    parser = argparse.ArgumentParser(
        prog="TLSM", description="""
       To read solar spectra on 1 nm resolution and output to the low-resolution
       binning scheme for NCAR/TGCM
       """
    )


    parser.add_argument(
        "-f","--bin_file",
        help="""Select binning scheme""",
        default=os.path.join(_rpath_,'input','low_bins.txt'),
        type = str
    )

    parser.add_argument(
        "-d","--indir",
        help="""Select input files directory""",
        default=os.path.join(_rpath_,'input','aia_masks', '2024', '099_new2'),
        type = str
    )

    parser.add_argument(
        "-p","--filepatt",
        help="""Select input matching pattern""",
        # default='fism_masked_*[0-9].nc*',
        default='FISM_MASKED_*[0-9].nc*',
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