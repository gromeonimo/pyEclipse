
import datetime
import netCDF4
import os
import sys
import numpy as np
import re


import ssl
import urllib
import urllib.request
import shutil


def ConvertDict(Config, verb):

    sections = Config.sections()
    rDict = {}

    for sect in sections:
        tDict = {}
        options = Config.options(sect)
        for option in options:
            try:
                tDict[option] = Config.get(sect, option)
                if tDict[option] == -1:
                    print("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                tDict[option] = None

        rDict[sect] = tDict

    return rDict

def get_spectra(solar_file, verb=0):
    """
 NAME:
 	get_spectra
 
 PURPOSE:
       Given an solar spectra file, read the solar spectra,
       convert wave length unit to Angstrom and solar flux 
       unit to photon cm^-2 s^-1 if necessary. 
 
 CATEGORY:
       utility function called by main program.
 
 CALLING SEQUENCE:
       out_spectra=get_spectra(solar_file)
 
 INPUTS:
       solar_file: (='dir_name/file_name')
 
 OUTPUTS:
       out_spectra: two dimensional array that holds 
                     output solar spectra. The array 
                     has at least 3 columns: waves and wavel 
                     in Angstrom, solar flux in
                     photon cm^-2 s^-1, any number of
                     extra flux columns. 
 
 PROCEDURE:
       some usual/useful conversion of units:
       1 nm = 10 A
       1 mW m^-2 = 1 erg cm^-2 s^-1
       flux(photon cm^-2 s^-1) = (6.242e11/12398)*wave(A) *flux(erg cm^-2 s^-1) 
 
 ROUTINES CALLED:
       read_gen, read_netcdf
 

    :param solar_file:   <input solar file>
    :param verb:         <verbosity level>
    :return:
    """
    if not os.path.exists(solar_file):
        print(f"ERROR - Unable to find file {solar_file}")
        return None

    pDict = {}
    pDict['netfile'] = solar_file
    pDict['ValuesOnly'] = 0
    pDict['getMetadata'] = 1
    pDict['reMapNames'] = 0
    jnk, nc_data = ReadNetcdfFile(pDict)

    if nc_data is None:
        print(f"ERROR - Unable to read input solar file {solar_file}")
        return None

    days = nc_data['variables']['DATE']['values'][0]
    n_days = len(days)


    w = nc_data['variables']['SP_WAVE']['values'][0]
    f = np.transpose(nc_data['variables']['SP_FLUX']['values'][0])
    n_lines = len(w)

    out_spectra = np.zeros((n_days + 1, n_lines), dtype=np.float64)
    out_spectra[0,:] = np.float64(w)
    out_spectra[1:, :] = np.float64(f.transpose())
    # for i in range(1,n_days+1):
    #     out_spectra[i,:] = np.float64(f[:,i-1])

    result = out_spectra.shape
    n_columns = result[0]

    # convert wave unit from nm to Angstrom, solar flux unit
    # from W M^-2 to photon cm^-2 s^-1 (from w m^-2 to erg cm^-2 s^-1 (1e3),
    # then from erg cm^-2 s^-1 to photon cm^-2 s^-1)

    out_spectra[0, :] = 10. * out_spectra[0, :]
    #out_spectra[1:, :] = (6.242e11 / 12398.) * out_spectra[0, :] * out_spectra[1:, :] * 1e3
    for i in range(1,n_columns):
        out_spectra[i, :] = (6.242e11 / 12398.0) * out_spectra[0, :]*out_spectra[i, :] * 1.0e3

    out_spectra = get_boundary(out_spectra, verb)

    return out_spectra


def get_spectra_fism2(solar_file, verb=0):
    """
 NAME:
 	get_spectra

 PURPOSE:
       Given an solar spectra file, read the solar spectra,
       convert wave length unit to Angstrom and solar flux
       unit to photon cm^-2 s^-1 if necessary.

 CATEGORY:
       utility function called by main program.

 CALLING SEQUENCE:
       out_spectra=get_spectra(solar_file)

 INPUTS:
       solar_file: (='dir_name/file_name')

 OUTPUTS:
       out_spectra: two dimensional array that holds
                     output solar spectra. The array
                     has at least 3 columns: waves and wavel
                     in Angstrom, solar flux in
                     photon cm^-2 s^-1, any number of
                     extra flux columns.

 PROCEDURE:
       some usual/useful conversion of units:
       1 nm = 10 A
       1 mW m^-2 = 1 erg cm^-2 s^-1
       flux(photon cm^-2 s^-1) = (6.242e11/12398)*wave(A) *flux(erg cm^-2 s^-1)

 ROUTINES CALLED:
       read_gen, read_netcdf


    :param solar_file:   <input solar file>
    :param verb:         <verbosity level>
    :return:
    """
    if not os.path.exists(solar_file):
        print(f"ERROR - Unable to find file {solar_file}")
        return None

    pDict = {}
    pDict['netfile'] = solar_file
    pDict['ValuesOnly'] = 0
    pDict['getMetadata'] = 1
    pDict['reMapNames'] = 0
    jnk, nc_data = ReadNetcdfFile(pDict)

    if nc_data is None:
        print(f"ERROR - Unable to read input solar file {solar_file}")
        return None

    #days = nc_data['variables']['DATE']['values'][0]
    days = np.asarray([int(nc_data['attributes']['Year_DOY'])])
    secs = np.asarray([int(nc_data['attributes']['SecondOfDay'])])
    n_days = len(days)
    date = None
    if n_days == 1:
        date = datetime.datetime.strptime(str(days[0]),"%Y%j")
        date += datetime.timedelta(seconds=int(secs[0]))

    w = nc_data['variables']['Wavelengths']['values']
    f = np.transpose(nc_data['variables']['msk_irradiance']['values'])
    n_lines = len(w)

    shape = list(f.shape)
    nlalo = shape[1] * shape[2]
    # newf = np.sum(f.reshape([shape[0], nlalo]), axis=1).reshape([shape[0], n_days])
    newf = f.reshape([shape[0], nlalo], order='C')
    n_days = nlalo


    out_spectra = np.zeros((n_days + 1, n_lines), dtype=np.float64)
    out_spectra[0, :] = np.float64(w)
    out_spectra[1:,:] = np.float64(newf[:, :].transpose())
    # for i in range(1, n_days + 1):
    #     out_spectra[i, :] = np.float64(newf[:, i - 1])

    result = out_spectra.shape
    n_columns = result[0]

    # convert wave unit from nm to Angstrom, solar flux unit
    # from W M^-2 to photon cm^-2 s^-1 (from w m^-2 to erg cm^-2 s^-1 (1e3),
    # then from erg cm^-2 s^-1 to photon cm^-2 s^-1)

    # out_spectra[1:, :] = (6.242e11 / 12398.) * out_spectra[0, :] * out_spectra[1:, :] #* 1e3
    out_spectra[1:, :] = out_spectra[1:, :] * 0.1 #* 0.175
    # for i in range(1, n_columns):
    #     out_spectra[i, :] = (6.242e11 / 12398.0) * out_spectra[0, :] * out_spectra[i, :] * 10. * 1.0e3 * 0.1




    out_spectra[0, :] = 10. * out_spectra[0, :]

    out_spectra = get_boundary(out_spectra, verb)

    data0 = out_spectra[0,:]
    idx0 = np.where(data0 < 0)[0]
    if len(idx0)> 0:
        data0[idx0] = 1.0e-10
        out_spectra[0, :] = data0



    return out_spectra, nc_data


def get_spectra_saveset(indata, verb=0):
    """
 NAME:
 	get_spectra

 PURPOSE:
       Given an solar spectra file, read the solar spectra,
       convert wave length unit to Angstrom and solar flux
       unit to photon cm^-2 s^-1 if necessary.

 CATEGORY:
       utility function called by main program.

 CALLING SEQUENCE:
       out_spectra=get_spectra(solar_file)

 INPUTS:
       solar_file: (='dir_name/file_name')

 OUTPUTS:
       out_spectra: two dimensional array that holds
                     output solar spectra. The array
                     has at least 3 columns: waves and wavel
                     in Angstrom, solar flux in
                     photon cm^-2 s^-1, any number of
                     extra flux columns.

 PROCEDURE:
       some usual/useful conversion of units:
       1 nm = 10 A
       1 mW m^-2 = 1 erg cm^-2 s^-1
       flux(photon cm^-2 s^-1) = (6.242e11/12398)*wave(A) *flux(erg cm^-2 s^-1)

 ROUTINES CALLED:
       read_gen, read_netcdf


    :param solar_file:   <input solar file>
    :param verb:         <verbosity level>
    :return:
    """


    if indata is None:
        print(f"ERROR - No data to process")
        return None


    irradiance = indata['irradiance']
    irr_shape = list(irradiance.shape)
    nmin = irr_shape[0]


    w = indata['wavelength']
    f = np.transpose(indata['irradiance'])
    n_lines = len(w)
    shape = list(f.shape)
    nlalo = shape[1]
    # newf = np.sum(f.reshape([shape[0], nlalo]), axis=1).reshape([shape[0], n_days])
    newf = f.reshape([shape[0], nlalo], order='C')
    n_days = nlalo

    out_spectra = np.zeros((n_days + 1, n_lines), dtype=np.float64)
    out_spectra[0, :] = np.float64(w)
    out_spectra[1:,:] = np.float64(newf[:, :].transpose())
    # for i in range(1, n_days + 1):
    #     out_spectra[i, :] = np.float64(newf[:, i - 1])

    result = out_spectra.shape
    n_columns = result[0]


    # convert wave unit from nm to Angstrom, solar flux unit
    # from W M^-2 to photon cm^-2 s^-1 (from w m^-2 to erg cm^-2 s^-1 (1e3),
    # then from erg cm^-2 s^-1 to photon cm^-2 s^-1)

    # out_spectra[1:, :] = (6.242e11 / 12398.) * out_spectra[0, :] * out_spectra[1:, :] #* 1e3
    # for i in range(1, n_columns):
    #     out_spectra[i, :] = (6.242e11 / 12398.0) * out_spectra[0, :] * out_spectra[i, :] * 10. * 1.0e3 * 0.1


    out_spectra[0, :] = 10. * out_spectra[0, :]
    out_spectra[1:, :] = out_spectra[1:, :] / 10.0

    out_spectra = get_boundary(out_spectra, verb)

    data0 = out_spectra[0,:]
    idx0 = np.where(data0 < 0)[0]
    if len(idx0)> 0:
        data0[idx0] = 1.0e-10
        out_spectra[0, :] = data0



    return out_spectra

    # #days = nc_data['variables']['DATE']['values'][0]
    # days = np.asarray([int(nc_data['attributes']['Year_DOY'])])
    # secs = np.asarray([int(nc_data['attributes']['SecondOfDay'])])
    # n_days = len(days)
    # date = None
    # if n_days == 1:
    #     date = datetime.datetime.strptime(str(days[0]),"%Y%j")
    #     date += datetime.timedelta(seconds=int(secs[0]))
    #
    # w = nc_data['variables']['Wavelengths']['values']
    # f = np.transpose(nc_data['variables']['msk_irradiance']['values'])
    # n_lines = len(w)
    #
    # shape = list(f.shape)
    # nlalo = shape[1] * shape[2]
    # # newf = np.sum(f.reshape([shape[0], nlalo]), axis=1).reshape([shape[0], n_days])
    # newf = f.reshape([shape[0], nlalo], order='C')
    # n_days = nlalo
    #
    #
    # out_spectra = np.zeros((n_days + 1, n_lines), dtype=np.float64)
    # out_spectra[0, :] = np.float64(w)
    # out_spectra[1:,:] = np.float64(newf[:, :].transpose())
    # # for i in range(1, n_days + 1):
    # #     out_spectra[i, :] = np.float64(newf[:, i - 1])
    #
    # result = out_spectra.shape
    # n_columns = result[0]
    #
    # # convert wave unit from nm to Angstrom, solar flux unit
    # # from W M^-2 to photon cm^-2 s^-1 (from w m^-2 to erg cm^-2 s^-1 (1e3),
    # # then from erg cm^-2 s^-1 to photon cm^-2 s^-1)
    #
    # # out_spectra[1:, :] = (6.242e11 / 12398.) * out_spectra[0, :] * out_spectra[1:, :] #* 1e3
    # # for i in range(1, n_columns):
    # #     out_spectra[i, :] = (6.242e11 / 12398.0) * out_spectra[0, :] * out_spectra[i, :] * 10. * 1.0e3 * 0.1
    #
    #
    # out_spectra[0, :] = 10. * out_spectra[0, :]
    #
    # out_spectra = get_boundary(out_spectra, verb)
    #
    # data0 = out_spectra[0,:]
    # idx0 = np.where(data0 < 0)[0]
    # if len(idx0)> 0:
    #     data0[idx0] = 1.0e-10
    #     out_spectra[0, :] = data0
    #
    #
    #
    # return out_spectra, nc_data




def get_boundary(data_in, verb=0):
    """
 +
  NAME:
        get_boundary
 
  PURPOSE:
        To calculate short boundary and long boundary 
        wave length for a centered wave length input spectra.
        Caution: for a non-uniformly spaced spectrum, it can cause
                 inconsistent result during the process of rebin.
 
  CATEGORY:
        utility function.
 
  CALLING SEQUENCE:
        data_out=get_boundary(data_in)
 
  INPUTS:
        data_in: two dimensional array that holds the input 
                 array. It has at least two columns:
                 (wave length, solar flux, scaling 
                 fluxes/additional fluxes if any)
      
  OUTPUTS:
        data_out: two dimensional array that holds output 
                  array. It has at least three columns: 
                  (short boundary,long boundary, solar flux, 
                  scaling fluxes/additional fluxes if any)
 
  COMMON BLOCKS:
        None.
 
  PROCEDURE:
 
  ROUTINES CALLED:
        None.
 
  MODIFICATION HISTORY:
        12/12/02, Liying Qian, Initial Version
 


    
    :param data_in:
    :param verb:
    :return:
    """
    wave = data_in[0,:]
    n_rows = len(wave)
    result = data_in.shape
    n_columns = result[0]
    data_out = np.zeros((n_columns+1,n_rows), dtype=np.float64)
    delw = np.zeros((n_rows-1), dtype=np.float64)
    waves = np.zeros(n_rows, dtype=np.float64)
    wavel = np.zeros(n_rows, dtype=np.float64)

    for i in range(n_rows-1):
        delw[i] = (wave[i+1]-wave[i])/2
        wavel[i] = wave[i]+delw[i]
        waves[i+1] = wavel[i]


    waves[0] = wave[0]-delw[0]
    wavel[n_rows-1] = wave[n_rows-1]+delw[n_rows-2]

    data_out[0,:] = waves
    data_out[1,:] = wavel
    data_out[2:,:] = data_in[1:,:]


    return data_out

def ReadNetcdfFile(pDict):

    infile = pDict['netfile']
    omsg = None
    if 'omsg' in pDict:
        omsg = pDict['omsg']
    verb = 0
    if verb in pDict:
        verb = pDict['verb']
    dtype = None
    if 'dtype' in pDict:
        dtype = pDict['dtype']

    getMetaData = 0
    if 'getMetadata' in pDict:
        getMetaData = pDict['getMetadata']

    reMapNames = 1
    if 'reMapNames' in pDict:
        reMapNames = pDict['reMapNames']

    if verb >= 3:
        if omsg is None:
            print('INFO - Parsing netCDF file '+infile)
        else:
            omsg.printout('Parsing netCDF file '+infile, 'info')

    if os.path.exists(infile):
        fname = os.path.basename(infile)
        try:
            ncd = netCDF4.Dataset(infile, 'r')
            ncd.set_auto_mask(False)
            values = {}
            attributes = {}
            for atname in ncd.ncattrs():
                attributes[atname] = ncd.getncattr(atname)

            dimensions = {}
            for dname in ncd.dimensions.keys():
                dimensions[dname] = ncd.dimensions[dname].size

            variables = {}

            for varname, variable in ncd.variables.items():
                #print("Dealing with "+fname+':'+varname)
                vars = {}
                vatts = {}
                for attrname in variable.ncattrs():
                    vatts[attrname] = variable.getncattr(attrname)
                vars['attributes'] = vatts
                dims = variable.dimensions
                ndims = len(dims)
                vars['dimensions'] = dims
                try:
                    Dtype = np.float64
                    tmptype = str(ncd.variables[varname].dtype)
                    if re.match(r'int', tmptype):
                        Dtype = np.int32
                    elif re.match(r'float32', tmptype):
                        Dtype = np.float32
                    else:
                        Dtype = None
                    if Dtype is not None:
                        vars['values'] = np.asarray(ncd.variables[varname][:], dtype=Dtype)
                    else:
                        vars['values'] = np.asarray(ncd.variables[varname][:])
                except:
                    if omsg is not None:
                        omsg.printout("Error reading variable {}".format(varname), 'error')
                    else:
                        print("Error reading variable {}".format(varname), file=sys.stderr)


                variables[varname] = vars


            values['attributes'] = attributes
            values['dimensions'] = dimensions
            values['variables']  = variables
            Groups = {}
            for g in ncd.groups:
                if g == 'Metadata' and getMetaData == 0:
                    continue

                tmp_grp = ncd.groups[g]
                gDict = {}

                aDict = {}
                for attr in tmp_grp.ncattrs():
                    aDict[attr] = tmp_grp.getncattr(attr)
                gDict['gattributes'] = aDict

                vDict = {}
                for varname, ncvar in ncd.groups[g].variables.items():
                    tDict = {}
                    if isinstance(ncvar[:], str):
                        tDict['values'] = ncvar[:]
                        tDict['dimensions'] = ncvar.dimensions
                        tDict['dtype'] = ncvar.dtype
                    else:
                        tDict['values'] = np.asarray(ncvar)
                        tDict['dimensions'] = ncvar.dimensions
                        tDict['dtype'] = ncvar.dtype
                    for attname in ncvar.ncattrs():
                        tDict[attname] = getattr(ncvar, attname)
                    vDict[varname] = tDict
                gDict['variables'] = vDict

                dDict = {}
                for dimname, dim in ncd.groups[g].dimensions.items():
                    dDict[dimname] = dim
                gDict['dimensions'] = dDict
                Groups[g] = gDict
            ng = len(Groups)
            if ng > 0:
                values['Groups'] = Groups




            ncd.close()
            rvalues = None
            # if reMapNames == 1:
            #     rvalues = reMapItems(values, infile, verb)

            if 'ValuesOnly' in pDict and pDict['ValuesOnly'] > 0:
                rDict = {}
                variables = rvalues['variables']
                for var in variables:
                    rDict[var] = variables[var]['values']

                return rDict, values


            return rvalues, values


        except:
            if omsg is None:
                print('ERROR - Unable to read file ' + infile)
            else:
                omsg.printout('Unable to read file  '+infile, 'error')

            return None, None
    else:
        if omsg is None:
            print('ERROR - Unable to find file ' + infile)
        else:
            omsg.printout('Unable to find file  '+infile, 'fatal')


def getbins(bin_file, verb):
    """
    Reads a binning file

    :param bin_file:
    :param verb:
    :return:
    """
    if verb > 1:
        print(f"INFO - Reading file {bin_file}")

    inf = None
    try:
        inf = open(bin_file, 'r')
    except:
        print(f"ERROR - Unable to tread bin file {bin_file}")
        sys.exit(1)
    lines = inf.readlines()
    inf.close()

    items = lines[1].split()
    nbins = int(items[0])
    scale = float(items[1])

    wave1 = np.zeros(nbins, dtype=np.float32)
    wave2 = np.zeros(nbins, dtype=np.float32)

    for i in range(2,nbins+2):
        items = lines[i].split()
        j = i - 2
        wave1[j] = float(items[0])
        wave2[j] = float(items[1])

    wave1 = wave1 * scale
    wave2 = wave2 * scale

    return wave1, wave2







def hi_to_nm(data_in, bpath, verb=0):
    """
    get 1 nm bins:

    :param data_in:
    :param verb:
    :return:
    """
    bin_file = os.path.join(bpath, 'input', 'nm_bins.txt')

    if not os.path.exists(bin_file):
        print(f"ERROR Unable to find bin file {bin_file}")
        return None



    wave1, wave2 = getbins(bin_file, verb)
    n_bins = len(wave1)

    waves = data_in[0,:]
    wavel = data_in[1,:]

    # construct output spectra data structure
    result = data_in.shape
    n_columns = result[0]
    n_lines = result[1]
    data_out = np.zeros((n_columns, n_bins), dtype=np.float64)
    data_out[0, :] = wave1
    data_out[1, :] = wave2
    n_cind = n_columns #- 1

    # map input spectra to 1 nm binning scheme
    # At this point, it is required that the input spectra
    # has at least three columns: wave short, wave long, solar flux,
    # additional columns (solar fluxes)

    count_ly = 0
    for i in range(0, n_bins):
        if (wave1[i] == 1215.67):  # special treatment for Ly-1
            index = np.where( (waves <= 1215.67) & (wavel >= 1215.67))[0]
            count_ly = len(index)
            if count_ly != 0:
                ly_1 = index[0]                  #  sometimes the original
                up_1 = index[0]-1                #  data bin that has ly_alpha
                down_1 = index[0]+1              #  solar flux might does not
                if up_1 >= 0:                    #  include 1215.67A, so search nearby
                    # modified to add if statement because TIMED/SEE version 9 include
	                #  negative data when data is missing
                    if data_in[2, ly_1] > 0:
                        temp = [data_in[2, ly_1], data_in[2, up_1]]
                    else:
                        temp = [data_in[20, ly_1], data_in[20, up_1]]

                    max_val = np.max(temp)
                    idx = np.where(temp == max_val)[0]
                    ii = idx[0]
                    if ii == 1:
                        ly_1 = up_1
                if down_1 <= (n_lines-1):
                    if data_in[2, ly_1] > 0:
                        temp = [data_in[2, ly_1], data_in[2, down_1]]
                    else:
                        temp = [data_in[20, ly_1], data_in[20, up_1]]

                    max_val = np.max(temp)
                    idx = np.where(temp == max_val)[0]
                    ii = idx[0]
                    if ii == 1:
                        ly_1 = up_1
                index[0] = ly_1
                data_out[2:n_cind, i] = data_in[2:n_cind, index[0]]

            # index[0] = ly_1
            # data_out[2:n_cind, i] = data_in[2:n_cind, index[0]]

        else:
            ind = np.where( (wavel > wave1[i]) & (waves < wave2[i]) )[0]
            ni = len(ind)
            if ni <= 0:
                data_out[2:n_cind, i] = 0
            elif ni == 1:
                # I am using linear interpration if data is coarse
                # than the binning scheme of nm_bins.txt.ver1 which is very likely at wavelength <3.2nm
                delw_1 = wave2[i] - wave1[i]
                delw_2 = wavel[ind[0]] - waves[ind[0]]
                if delw_1 <= delw_2:
                    data_out[2:n_cind, i] = data_in[2:n_cind, ind[0]] * delw_1 / delw_2
                else:
                    data_out[2:n_cind, i] = data_in[2:n_cind, ind[0]]
            else:
                for j in range(0, ni):
                    # again, I am using linear interpolation
                    if ((waves[ind[j]] <= wave1[i]) and (wavel[ind[j]] > wave1[i])):
                        data_out[2:n_cind, i] = data_out[2:n_cind, i] + (wavel[ind[j]] - wave1[i]) / \
                                                (wavel[ind[j]] - waves[ind[j]]) * data_in[2:n_cind, ind[j]]
                    elif ( (waves[ind[j]] >= wave1[i]) and (wavel[ind[j]] <= wave2[i]) ):
                        data_out[2:n_cind, i] = data_out[2:n_cind, i] + data_in[2:n_cind, ind[j]]
                    elif ( (waves[ind[j]] < wave2[i]) and (wavel[ind[j]] >= wave2[i]) ):
                        data_out[2:n_cind, i] = data_out[2:n_cind, i] + (wave2[i] - waves[ind[j]]) / \
                                                (wavel[ind[j]] - waves[ind[j]]) * data_in[2:n_cind, ind[j]]

        if (wave1[i] == 1200 and count_ly != 0):
            data_out[2:n_cind, i] = data_out[2:n_cind, i] - data_in[2:n_cind, index[0]]



    if verb > 2:
        print(f"INFO - In function hi_to_nm, bin file is {bin_file}")
        print(data_out[0,:])
        print(data_out[1,:])
        print(data_out[50,:])

        # energy conservation check

        if n_columns == 3:
            print(f"in hi_to_nm, input: {np.sum(data_in[2, np.where(wavel <= 1750)[0]])}")
            print(f"in hi_to_nm, output: {np.sum(data_out[2, np.where(wave2 <= 1750)[0]])}")
        elif n_columns > 3:
            print(f"in hi_to_nm, input: {np.sum(data_in[2, np.where(wavel <= 1750)[0]]), np.sum(data_in[3,np.where(wavel <= 1750)[0]])}")
            print(f"in hi_to_nm, output: {np.sum(data_out[2,:]), np.sum(data_out[3,:])}")


            
    return data_out


def nm_to_hs(data_in, bpath, verb=0):
    """
  NAME:
        nm_to_hs
 
  PURPOSE:
        To transfer spectra on 1 nm binning scheme to Hinteregger reference
        spectrum resolution for wave length 100A-1000 A and
        add desired lines in FUV.
 
  CATEGORY:
        utility function called by the main program.
 
  CALLING SEQUENCE:
        data_out=nm_to_hs(data_in)
 
  INPUTS:
        data_in: a two dimension array that holds input
                 spectra on 1 nm bins. It has at least 3 
                 columns (short bounday in A, long boundary in A, 
                 solar flux in photon cm^-2 s^-1, additional
                 columns that are flux (photon cm^-2 s^-1) 
                 in nature if any).
 
  OUTPUTS:
        data_out: a two dimensional array that holds 
                  high resolution output spectra. It has same
                  data structure as data_in.
 
  PROCEDURE:
        The routine uses Hinteregger solar minimum reference 
        spectrum, calculates flux ratio of each wave length in
        the 1 nm it belongs to, then apply the ratio to the input 
        spectra to get spectra in Hinteregger resoltion (only
        for wave length 100A-1000A and desired lines in FUV). 
        This process conserves energy. 
 
  ROUTINES CALLED:
        read_hs  
 
    :param data_in:
    :param verb:
    :return:
    """
    # read hintergger reference spectrum
    hs_file = os.path.join(bpath, 'input', 'sc21refw.dat')

    if not os.path.exists(hs_file):
        print(f"ERROR Unable to find bin file {hs_file}")
        return None

    wave1 = data_in[0, :]
    wave2 = data_in[1, :]
    n_bins = len(wave1)
    result = data_in.shape
    n_columns = result[0]


    solar_activity='low'
    a = read_hs(hs_file,solar_activity,1661, verb)
    refwvln = a[0,:]
    refflux = a[2,:]
    n_hs = len(refwvln)

    #  the following arrays are used to hold the
    #  mapped data in the Hinteregger wave space
    #  Since Hinteregger spectrum does not have
    #  data for the following 4 nm bins (370-380A,
    #  380A-390A,420-430A,440-450A), the array are made 
    #  bigger than Hinteregger array
    n_cind = n_columns-2
    ratio = np.zeros(n_hs+10, dtype=np.float32)
    flux = np.zeros((n_cind,n_hs+10), dtype=np.float32)
    waves = np.zeros(n_hs+10, dtype=np.float32)
    wavel = np.zeros(n_hs+10, dtype=np.float32)

    #  Get ratio for each Hinteregger wave length.
    #  Since Hinteregger spectrum has poor quality 
    #  in FUV, here HS is only referenced for wave length
    #  50A-1050A, and it is enough for our calculation
    #  purpose for low resolution scheme, which has 
    #  logical bins in the wave length range 
    #  650-975A

    i_start = np.where( (wave1 <= 50) & (wave2 > 50))[0]
    scount = len(i_start)
    i_end = np.where( (wave1 <= 1050) & (wave2 > 1050))[0]
    ecount = len(i_end)
    if scount == 0: i_start[0]=0
    if ecount == 0: i_end[0]=n_bins-1

    index = 0     # number of Hinteregger spectrum that a ratio is calculated
                  # plus number of 1 nm bins that has no Hinteregger data

    for i in range(i_start[0], i_end[0]+1):
        ind = np.where((refwvln >= wave1[i]) & (refwvln < wave2[i]))[0]
        n_ind = len(ind)
        if n_ind > 0:
            tot_flux = np.sum(refflux[ind])
            if tot_flux == 0: tot_flux=0.001
            for k in range(0, n_ind):
                ratio[index] = refflux[ind[k]] / tot_flux
                #flux[0:n_cind - 1, index] = ratio[index] * data_in[2:n_columns - 1, i]
                flux[0:n_cind, index] = ratio[index] * data_in[2:n_columns, i]
                waves[index] = refwvln[ind[k]]
                wavel[index] = waves[index]
                index = index + 1
        else:
            waves[index] = data_in[0, i]
            wavel[index] = data_in[1, i]
            #flux[0:n_cind - 1, index] = data_in[2:n_columns - 1, i]
            flux[0:n_cind, index] = data_in[2:n_columns, i]
            index = index + 1

    #  the data_in is now mapped to Hinteregger wave length 
    #  space from i_start to i_end
    #  The output spectra will consists of its original spectra
    #  for wave length shortward of i_start and longward of
    #  i_end, and mapped spectra on Hinteregger space from
    #  i_start to i_end 
    data_out = None
    if i_end[0] == n_bins-1:
        n_out = i_start[0] + index
        data_out = np.zeros((n_columns, n_out), dtype=np.float32)
    else:
        n_out = i_start[0] + index + (n_bins - 1 - i_end[0])
        data_out = np.zeros((n_columns, n_out), dtype=np.float32)
        data_out[:, i_start[0] + index:n_out] = data_in[:, i_end[0] + 1:n_bins]

    if i_start[0] > 0: data_out[:, 0:i_start[0]] = data_in[:, 0:i_start[0]]
    for i in range(0, index):
        data_out[0, i_start[0] + i] = waves[i]
        data_out[1, i_start[0] + i] = wavel[i]
        data_out[2:n_columns, i_start[0] + i] = flux[0:n_cind, i]

    #  GLOW and low resolution scheme have lines. Lines less than 1000A
    #  have been taken care in the above step. Lines in FUV need to add 
    #  to the spectra as well.

    #  energy conservation check

    if verb > 1:
        if n_columns == 3:
            print(f"in nm_to_hs, input:  {np.sum(data_in[2, :])}")
            print(f"in nm_to_hs, output: {np.sum(data_out[2, :])}")
        elif n_columns > 3:
            print(f"in nm_to_hs, input:  {np.sum(data_in[2, :])}, {np.sum(data_in[3, :])}")
            print(f"in nm_to_hs, output: {np.sum(data_out[2,:])}, {np.sum(data_out[3,:])}")

        print()


    return data_out




def read_hs(solar_file,solar_activity,n_rows, verb=0):
    """
  NAME:
        read_hs
 
  PURPOSE:
        to read Hinteregger Spectrum.
 
  CATEGORY:
        utility function.
 
  CALLING SEQUENCE:
        data=read_hs(solar_file,solar_activity,n_rows)
 
  INPUTS:
        solar_file:     source file name of solar spectrum.
        solar_activity: 'low' for solar minimum,
                        'high' for solar maximum.
        n_rows:          number of lines of data.
 
  OUTPUTS:
        data:  a two dimensional array that holds spectrum data.
               It has 4 columns for solar minimum (wave length,
               solar flux, origin of solar irradiance, flux ratio).
               It has 2 columns for solar maximum (wave length,
               solar flux).
 
  COMMON BLOCKS:
        None.
  
  PROCEDURE:
        This routine modifies the original solar spectrum to include 
        a line (789.36 A) in the 5 nm bins and lines' binning scheme.
     
  ROUTINES CALLED:
        None.
 
  MODIFICATION HISTORY:
        12/12/02, Liying Qian, Initial version.


    :param solar_file:
    :param solar_activity:
    :param n_rows:
    :param verb:
    :return:
    """

    temp = np.zeros((4,n_rows), dtype=np.float32)

    if verb > 1:
        print(f"INFO - Reading file {solar_file}")

    inf = None
    try:
        inf = open(solar_file, 'r')
    except:
        print(f"ERROR - Unable to tread bin file {solar_file}")
        sys.exit(1)
    lines = inf.readlines()
    inf.close()

    for i in range(0, n_rows):
        items = lines[i].split()
        if solar_activity == 'low':
            temp[0,i] = float(items[0])
            temp[1,i] = float(items[1]) * 1.0e6  # to be in units of photons/cm2/s
            temp[2,i] = int(items[-2])
            temp[3,i] = float(items[-1])
        else:
            temp[0,i] = float(items[0])
            temp[1,i] = float(items[1]) * 1.0e6  # to be in units of photons/cm2/s

    # modify the original solar spectrum to include a line in GLOW at
    # 789.36A

    ind = np.where((temp[0, :] == 787.71) | (temp[0, :] == 790.15))[0]
    temp_flux = np.sum(temp[1, ind])
    temp[1, ind] = 0.0
    temp_ratio = temp[3, ind[0]]
    temp[3, ind] = 0.0
    n_temp = len(temp[0, :])
    b = np.zeros((4, n_temp + 1), dtype=np.float32)
    b[:, 0: n_temp] = temp
    b[0, n_temp] = 789.36
    b[1, n_temp] = temp_flux
    b[3, n_temp] = temp_ratio
    b[2, n_temp] = temp[2, ind[0]]
    ind = np.argsort(b[0, :])
    b = b[:, ind]

    # Assign data for the short wave length where data is not
    # available.

    # data at short wave length
    arr_l = np.asarray([[1, 1e-1, 2, 1e-1], [2, 5e1, 2, 5e1], [4, 1e4, 2, 1e4], [8, 2e6, 2, 2e6]])
    arr_h = np.asarray([[1, 5e2, 0, 0], [2, 3e4, 0, 0], [4, 8e5, 0, 0], [8, 5e7, 0, 0]])

    # append it to the original data

    # arr_l = np.transpose(arr_l)
    # arr_h = np.transpose(arr_h)
    b = np.transpose(b)
    b_org = np.copy(b)
    if solar_activity == 'low':
        #b = np.concatenate((arr_l, b), axis=0)
        b = np.vstack((arr_l, b))
    else:
        b = np.concatenate((arr_h, b), axis=0)

    b = np.transpose(b)
    n_rows = len(b[0, :])

    # convert scaling indices to flux to adapt the data to 
    # general procedures of TLSM

    if solar_activity == 'low':
        sc1 = np.zeros(n_rows, dtype=np.float32)
        sc2 = np.zeros(n_rows, dtype=np.float32)

        for i in range(n_rows):
            if b[2, i] == 1:
                sc2[i] = b[3, i] * b[1, i]
                sc1[i] = 0
            elif b[2,i] == 2:
                sc2[i] = b[3, i] * b[1, i]
                sc1[i] = 0

        b[2, :] = sc1
        b[3, :] = sc2

    # change the format to adapt the spectrum to general 
    # procedure rebin
    data = np.zeros((5, n_rows), dtype=np.float32)
    data[0, :] = b[0, :]
    data[1, :] = b[0, :]
    data[1, 0:4] = [2., 4., 8., 16.]
    data[2:5, :] = b[1:4, :]

    # the scale2 value for the first 4 bins (1-16A) are obtained from GLOW.
    # in GLOW, scale1*1.e6 would be in units of photon cm^-2 s^-1, here in
    # this file, scale1 and scale2 should be in units of photon cm^-2 s^-1

    data[4, 0] = 4.8E-6 * 1.e6
    data[4, 1] = 2.9E-4 * 1.e6
    data[4, 2] = 7.6E-3 * 1.e6
    data[4, 3] = 0.46 * 1.e6

    return data

def rebin(bin_file, in_spectra, indir, verb=0):
    """
 NAME:
       rebin

 PURPOSE:
       To put input solar spectra to a selected binning scheme

 CATEGORY:
       core procedure called by the main program.

 CALLING SEQUENCE:
       rebin,bin_file,in_spectra,out_spectra,sigma_o,sigma_o2,sigma_n,sigma_n2

 INPUTS:
       bin_file:     the selected binning scheme file.
       in_spectra:   a two dimensional array that holds the input
                     spectra. It has at least 3 columns:
                     wave short boundary in A, long boundary in A, 
                     solar flux in photon cm^-2 s^-1, more fluxes if any.
                     For example, for Hinteregger solar minimum
                     spectrum, in_spectra has 5 columns:
                     wave short, wave long, solar flux, flux from
                     chromsphere, flux from corona. For Woods 
                     reference spectrum, it also has 5 columns:
                     wave short, wave long, solar flux, 
                     27-day variation, 11-year variation.

 OUTPUTS:
       out_spectra:  a two dimensional array that holds the output
                     spectra. It has same data structure as in_spectra.

 ROUTINES CALLED:
       getbins, read_gen



    :param bin_file:
    :param in_spectra:
    :param verb:
    :return:
    """

    wave1, wave2 = getbins(bin_file, verb)
    n_bins = len(wave1)

    iflines = (wave1 == wave2) * 1    # 1 where bins are lines, 0 for continuum bins

    # unpack input spectra

    refwvln = 0.5 * (in_spectra[0, :]+in_spectra[1, :])
    n_rows = len(refwvln)
    result = in_spectra.shape
    n_columns = result[0]
    n_cind = n_columns - 2
    refflux = np.zeros((n_cind, n_rows), dtype=np.float32)
    refflux[0:n_cind, :] = in_spectra[2:n_columns, :]

    # check lines in the input spectra
    lines = wave1 * iflines
    ifreflines = np.zeros(n_rows, dtype=np.float32)
    for i in range(0, n_rows):
        for j in range(0,n_bins):
            if ((lines[j] == refwvln[i]) and (iflines[j] == 1)):
                ifreflines[i] = 1

    # read cross section data. Logical bins are divided based on N2 absorption 
    # cross sections

    fenn = read_gen(os.path.join(indir, 'input', 'phfenn.tab'), 1, 10, 1945, verb)
    henken = read_gen(os.path.join(indir, 'input', 'henken.tab'), 1, 3, 342, verb)

    # Combine Fennely, Henke
    ind_h2 = np.where(henken[0, :] <= 50)[0]
    ind_f2 = np.where(fenn[0, :] > 50)[0]
    a1 = None
    b1 = None
    a2 = None
    b2 = None
    if len(ind_h2)> 0:
        a1 = np.transpose(henken[0,ind_h2])
        a2 = np.transpose(henken[2,ind_h2])
    if len(ind_f2) > 0:
        b1 = np.transpose(fenn[0,ind_f2])
        b2 = np.transpose(fenn[1,ind_f2])

    if a1 is not None and b1 is not None:
        ln2 = np.concatenate([a1, b1])
    elif a1 is None and b1 is not None:
        ln2 = np.copy(b1)
    elif a1 is not None and b1 is None:
        ln2 = np.copy(a1)


    if a2 is not None and b2 is not None:
        an2_ext = np.concatenate([a2, b2])
    elif a2 is None and b2 is not None:
        an2_ext = np.copy(b2)
    elif a2 is not None and b2 is None:
        an2_ext = np.copy(a2)

    fennsigan2 = np.interp(refwvln, ln2, an2_ext)

    # Get solar flux into model bins.

    modflux = np.zeros((n_cind, n_bins), dtype=np.float32)

    flag_low = (n_bins < 45)
    count_1 = 0
    count_2 = 0
    count_3 = 0

    for ibin in range(0, n_bins):
        # the following if statement should handle lines or continua, it does
        # assume though that the line is one of the Hinteregger lines

        if iflines[ibin] == 1:
            ind = np.where((refwvln >= wave1[ibin]) & (refwvln <= wave2[ibin]))[0]
        else:
            if flag_low == 1:
                match wave1[ibin]:
                    case 650.0:
                        if count_1 == 0:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin]) & \
                                           (fennsigan2 < 31) )[0]
                            count_1 = count_1 + 1
                        elif count_1 == 1:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin]) & \
                                           (fennsigan2 >= 31) )[0]

                    case 798.0:
                        if count_2 == 0:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin]) &  \
                                           (fennsigan2 < 4) )[0]
                            count_2 = count_2 + 1
                        elif count_2 == 1:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin]) & \
                                           (fennsigan2 >= 4) & (fennsigan2 < 31) )[0]
                            count_2 = count_2 + 1
                        elif count_2 == 2:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin])  & \
                                           (fennsigan2 >= 31) )[0]

                    case 913.0:
                        if count_3 == 0:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin])  & \
                                           (fennsigan2 < 4) )[0]
                            count_3 = count_3 + 1
                        elif count_3 == 1:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin])  & \
                                           (fennsigan2 >= 4) & (fennsigan2 < 31) )[0]
                            count_3 = count_3 + 1
                        elif count_3 == 2:
                            ind = np.where((refwvln >= wave1[ibin]) & \
                                           (refwvln < wave2[ibin])  & \
                                           (fennsigan2 >= 31) )[0]

                    case _:
                        ind = np.where((refwvln >= wave1[ibin]) & (refwvln < wave2[ibin]))[0]

            else:
                ind = np.where((refwvln >= wave1[ibin]) & (refwvln < wave2[ibin]))[0]

            ni = len(ind)
            if ni > 0:
                indcont = np.where(ifreflines[ind] == 0)[0]
                oldind = np.copy(ind)
                nid = len(indcont)
                if nid > 0:
                    indcont = ind[indcont]
                    ind = indcont

        ni = len(ind)
        if ni > 0:
            for cind in range(0, n_cind):
                modflux[cind, ibin] = np.sum(refflux[cind, ind])

    # output

    out_spectra = np.zeros((n_columns, n_bins), dtype=np.float32)
    out_spectra[0, :] = wave1
    out_spectra[1, :] = wave2
    out_spectra[2:n_columns, :] = modflux[0:n_cind, :]

    # energy conservation check
    if verb > 1:
        ind1 = np.where(in_spectra[0, :] < wave2[n_bins - 1])[0]
        ind2 = np.where(refwvln < wave2[n_bins - 1])[0]

        if n_columns == 3:
            print(f"in rebin, input:  {np.sum(in_spectra[2, ind1])}")
            print(f"in rebin, output: {np.sum(out_spectra[2, :])}")
        elif n_columns > 3:
            print(f"in rebin, input:  {np.sum(in_spectra[2, ind1])}, {np.sum(in_spectra[3, ind1])}")
            print(f"in rebin, output: {np.sum(out_spectra[2, :])}, {np.sum(out_spectra[3, :])}")

    return out_spectra




def read_gen(filename,n_comm,n_columns,n_lines, verb=0):
    if not os.path.exists(filename):
        print(f"ERROR - In read-gen, unable to find file {filename}")
        sys.exit(1)

    data = np.zeros((n_columns, n_lines), dtype=np.float32)
    inf = open(filename, 'r')
    lines = inf.readlines()
    inf.close()
    for i in range(n_comm, n_lines+1):
        items = lines[i].split()
        data[:,i-n_comm] = [float(it) for it in items]

    return data

def put_spectra(solar_file, out_spectra, verb=0):
    """
 NAME:
	put_spectra

 PURPOSE:
       To write outputs

 CATEGORY:
       utility procedure called by main program.

 CALLING SEQUENCE:
       put_spectra,solar_file,out_spectra

 INPUTS:
       solar_file:  (='dir_name/file_name'), original input
                    spectra file name
       out_spectra: output spectra 


 PROCEDURE:

 ROUTINES CALLED:
       None.

    :param solar_file:
    :param out_spectra:
    :param verb:
    :return:
    """

    shape = out_spectra.shape
    n_bins = shape[1]

    sfname = os.path.basename(solar_file)
    indir = os.path.dirname(solar_file)
    if re.search(r'\.ncdf$', sfname, re.IGNORECASE) or re.search(r'\.nc4{0,1}', sfname, re.IGNORECASE):

        n_days = shape[0] - 2
        pDict = {}
        pDict['netfile'] = solar_file
        pDict['ValuesOnly'] = 0
        pDict['getMetadata'] = 1
        pDict['reMapNames'] = 0
        jnk, data_old = ReadNetcdfFile(pDict)

        data_new = {}
        UTTIME = np.zeros(n_days, dtype=np.float64)
        if re.search(r'L3A', sfname, re.IGNORECASE) or re.search(r'FISM', sfname, re.IGNORECASE):
            if 'TIME' in data_old['variables']:
                UTTIME = data_old['variables']['TIME']['values']/3600.
            elif 'SecondOfDay' in data_old['attributes']:
                UTTIME = data_old['attributes']['SecondOfDay']/3600.

        else:
            UTTIME[:] = 12.0

        att_file = os.path.join(indir, 'see__L3.attr')
        if not os.path.exists(att_file):
            att_file = None

        data_new['UTTIME'] = UTTIME

        yfrac = np.zeros(n_days, dtype=np.float64)
        dates = data_old['variables']['DATE']['values'][0]
        for i in range(0, n_days):
            iyear = int(dates[i] / 1000)
            if iyear % 4  == 0:
                yfrac[i] = int(dates[i] / 1000) + (float(dates[i]-iyear*1000)-1 + float(UTTIME[i]/24.))/366.
            else:
                yfrac[i] = int(dates[i] / 1000) + (float(dates[i]-iyear*1000)-1 + float(UTTIME[i]/24.))/365.

        data_new['YFRAC'] = yfrac
        data_new['DATE'] = dates


        temp_array=np.flip(out_spectra,1)
        data_new['WAVE1'] = temp_array[0,:]
        data_new['WAVE2'] = temp_array[1,:]
        data_new['SP_FLUX'] = temp_array[2:, :]

        dims = {'ndays':n_days, 'nbins':n_bins}

        name = sfname.split('.')[0]
        outfile = os.path.join(indir, name+'_'+str(n_bins)+'_py.nc')
        if os.path.exists(outfile):
            os.remove(outfile)


        if verb:
            print(f"INFO - Writing to file {outfile}")

        write_netcdf(data_new, outfile, dims, att_file, verb)



def put_spectra_fism2(solar_file, out_spectra, init_data, verb=0):
    """
 NAME:
	put_spectra

 PURPOSE:
       To write outputs

 CATEGORY:
       utility procedure called by main program.

 CALLING SEQUENCE:
       put_spectra_fism2,solar_file,out_spectra,init_data

 INPUTS:
       solar_file:  (='dir_name/file_name'), original input
                    spectra file name
       out_spectra: output spectra


 PROCEDURE:

 ROUTINES CALLED:
       None.

    :param solar_file:
    :param out_spectra:
    :param verb:
    :return:
    """
    nLat = None
    nLon = None
    try:
        nLat = init_data['dimensions']['nLat']
        nLon = init_data['dimensions']['nLon']
    except:
        print(f"ERROR - Unable to get Lat/Lon dimensions from file {solar_file}")
        return


    shape = out_spectra.shape
    n_bins = shape[1]

    start_wave = out_spectra[0,:]#/10.0
    end_wave = out_spectra[1,:]#/10.0
    # msk_irradiance = (out_spectra[2:, :].transpose()).reshape([n_bins, nLat, nLon])
    #msk_irradiance = (out_spectra[2:, :]).reshape([nLat, nLon, n_bins], order='C')
    msk_irradiance = (out_spectra[2:, :]).reshape([nLat, nLon, n_bins], order='F')

    #tempdata = (out_spectra[2:, :].transpose()).reshape([n_bins, nLat, nLon], order='C')



    # wavelengths = (start_wave + end_wave)/2.0
    # deltawave = (end_wave - start_wave)/2.0
    wavelengths = start_wave/10.0
    #deltawave = (end_wave - start_wave)/10.0
    #deltawave = (end_wave - start_wave)/10.0
    endwaves = end_wave/10.0




    sfname = os.path.basename(solar_file)
    indir = os.path.dirname(solar_file)

    items = sfname.split('.')
    ni = len(items)
    newname = ''
    for i in range(ni-1):
        newname += items[i]

    newname += '_TLSM.'+items[-1]

    outfile = os.path.join(indir, newname)
    if os.path.exists(outfile):
        if verb > 1:
            print(f"INFO - Removing file {outfile} to allow overwriting")
        os.remove(outfile)

    now = datetime.datetime.now()
    compress = 'zlib'
    ncdf = None
    try:
        ncdf = netCDF4.Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    except:
        print(f"ERROR - Unable to open file {outfile} for writing")
        return


    ncdf.CreationTime = now.strftime("%Y-%m-%d %H:%M:%S.%f")
    try:
        ncdf.Year_DOY = init_data['attributes']['Year_DOY']
        ncdf.SecondOfDay = init_data['attributes']['SecondOfDay']
    except:
        print(f"ERROR - Unable to add global attributes")


    dim_lat = ncdf.createDimension('nLat', nLat)
    dim_lon = ncdf.createDimension('nLon', nLon)
    dim_wav = ncdf.createDimension('nWave', n_bins)
    dim_wav2 = ncdf.createDimension('nWave2', init_data['dimensions']['nWave'])



    if 'iseclipse' in init_data['variables']:
        var = ncdf.createVariable('iseclipse', (np.float32), [dim_lat, dim_lon], zlib=1, complevel=5)
        var[:,:] = init_data['variables']['iseclipse']['values']

    var = ncdf.createVariable('Lat', (np.float32), dim_lat, zlib=1, complevel=5)
    var[:] = init_data['variables']['Lat']['values']

    var = ncdf.createVariable('Lon', (np.float32), dim_lon, zlib=1, complevel=5)
    var[:] = init_data['variables']['Lon']['values']

    var = ncdf.createVariable('msk_irradiance_tlsm', (np.float64), [dim_lat, dim_lon, dim_wav], zlib=1, complevel=5)
    var[:,:,:] = msk_irradiance[:,:,:]

    var = ncdf.createVariable('Start_wavelength', (np.float32), dim_wav, zlib=1, complevel=5)
    var[:] = wavelengths

    var = ncdf.createVariable('End_wavelength', (np.float32), dim_wav, zlib=1, complevel=5)
    var[:] = endwaves

    if 'sza' in init_data['variables']:
        var = ncdf.createVariable('sza', (np.float32), [dim_lat, dim_lon], zlib=1, complevel=5)
        var[:,:] = init_data['variables']['sza']['values']


    varfis = ncdf.createVariable('fism2', (np.float32), (dim_wav2), zlib=1, complevel=5)
    varfis[:] = init_data['variables']['fism2']['values']
    varwave = ncdf.createVariable('Wavelengths', (np.float32), (dim_wav2), zlib=1, complevel=5)
    varwave[:] = init_data['variables']['Wavelengths']['values']


    ncdf.close()



    # if verb:
    print(f"INFO - Wrote to file {outfile}")


def write_spectra_saveset(out_spectra, init_data, tag, verb=0):
    """
 NAME:
	put_spectra

 PURPOSE:
       To write outputs

 CATEGORY:
       utility procedure called by main program.

 CALLING SEQUENCE:
       put_spectra_fism2,solar_file,out_spectra,init_data

 INPUTS:
       solar_file:  (='dir_name/file_name'), original input
                    spectra file name
       out_spectra: output spectra


 PROCEDURE:

 ROUTINES CALLED:
       None.

    :param solar_file:
    :param out_spectra:
    :param verb:
    :return:
    """
    nLat = None
    nLon = None

    solar_file = init_data['savefile']

    shape = out_spectra.shape
    n_bins = shape[1]
    n_mins = shape[0] - 2

    start_wave = out_spectra[0,:]#/10.0
    end_wave = out_spectra[1,:]#/10.0
    # msk_irradiance = (out_spectra[2:, :].transpose()).reshape([n_bins, nLat, nLon])
    #msk_irradiance = (out_spectra[2:, :]).reshape([nLat, nLon, n_bins], order='C')
    msk_irradiance = (out_spectra[2:, :]).reshape([n_mins, n_bins], order='F')

    #tempdata = (out_spectra[2:, :].transpose()).reshape([n_bins, nLat, nLon], order='C')



    # wavelengths = (start_wave + end_wave)/2.0
    # deltawave = (end_wave - start_wave)/2.0
    wavelengths = start_wave
    endwaves = end_wave
    if tag == 'TLSM':
        wavelengths = start_wave/10.0
        endwaves = end_wave/10.0




    sfname = os.path.basename(solar_file)
    indir = os.path.dirname(solar_file)

    items = sfname.split('.')
    ni = len(items)
    newname = ''
    for i in range(ni-1):
        newname += items[i]

    newname += '_'+tag+'.nc4'

    outfile = os.path.join(indir, newname)
    if os.path.exists(outfile):
        if verb:
            print(f"INFO - Removing file {outfile} to allow overwriting")
        os.remove(outfile)

    now = datetime.datetime.now()
    compress = 'zlib'
    ncdf = None
    try:
        ncdf = netCDF4.Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
    except:
        print(f"ERROR - Unable to open file {outfile} for writing")
        return


    ncdf.CreationTime = now.strftime("%Y-%m-%d %H:%M:%S.%f")
    try:
        ncdf.Year_DOY = init_data['ydoy']
    except:
        print(f"ERROR - Unable to add global attributes")

    dim_min = ncdf.createDimension('n_minute', nLon)
    dim_wav = ncdf.createDimension('nWave', n_bins)


    var = ncdf.createVariable('irradiance', (np.float64), [dim_min, dim_wav], zlib=True)
    var[:,:] = msk_irradiance[:,:]

    var = ncdf.createVariable('Start_wavelength', (np.float32), dim_wav, zlib=True)
    var[:] = wavelengths

    var = ncdf.createVariable('End_wavelength', (np.float32), dim_wav, zlib=True)
    var[:] = endwaves

    var = ncdf.createVariable('utc', (np.float32), dim_min, zlib=True)
    var[:] = init_data['utc']

    var = ncdf.createVariable('jd', (np.float32), dim_min, zlib=True)
    var[:] = init_data['jd']

    ncdf.close()



    # if verb:
    print(f"INFO - Wrote to file {outfile}")



def write_netcdf(data, filename, dims, att_file=None, verb=0):
    """
 NAME:
	write_netCDF.pro

 PURPOSE:
	Write netCDF file given a structure variable

 CATEGORY:
	All levels of processing

 CALLING SEQUENCE:  
	write_netCDF, data, filename, status, path=dir_path, att_file=att_filename, /clobber

 INPUTS:
	data = structure variable of input data
	filename = filename for new netCDF file
	path = optional directory path for the attributes definition file
	att_file = optional filename for the attributes definition file
	clobber = optional option for creating netCDF file
			clobber means any old file will be destroyed

	An external *.att file is used to define attributes (where * = "data" structure name)

 OUTPUTS:  
	status = result status: 0 = OK_STATUS, -1 = BAD_PARAMS, -2 = BAD_FILE,
			-3 = BAD_FILE_DATA, -4 = FILE_ALREADY_OPENED

	A netCDF file is created and written.

 COMMON BLOCKS:
	None

 PROCEDURE:
	Check for valid input parameters
	Open the netCDF file
	Use the structure's tag names for defining the variable names in the netCDF.
	Use the structure name and optional 'path' variable for the Attributes filename
		OR use the optional 'att_file' parameter for this filename
	If this Attributes definition file exists, then transfer those attributes into the netCDF file
		OR else don't write any attributes to the netCDF file.
	Once netCDF variables and attributes are defined, then write the structure's data to netCDF file
	Close the netCDF file

	NetCDF IDL Procedures / Process:
	1. NCDF_CREATE: Call this procedure to begin creating a new file. The new file is put into define mode.
	2. NCDF_DIMDEF: Create dimensions for the file.
 	3. NCDF_VARDEF: Define the variables to be used in the file.
	4. NCDF_ATTPUT: Optionally, use attributes to describe the data.  Global attributes also allowed.
	4. NCDF_CONTROL, /ENDEF: Leave define mode and enter data mode.
	5. NCDF_VARPUT: Write the appropriate data to the netCDF file.
	6. NCDF_CLOSE: Close the file.

 MODIFICATION HISTORY:
	9/20/99		Tom Woods		Original release code, Version 1.00


    :param data_new:
    :param outfile:
    :param att_file:
    :param verb:
    :return:
    """

    # Generic "status" values
    OK_STATUS = 0
    BAD_PARAMS = -1
    BAD_FILE = -2
    BAD_FILE_DATA = -3
    FILE_ALREADY_OPENED = -4

    debug_mode = 0          # set to >= 1 if want to debug this procedure
                            # set to 2 if want to debug and force directory to special Woods Mac directory

    # check for valid parameters
    status = BAD_PARAMS
    if data is None and filename is None:
        print(f"ERROR - USAGE: write_netCDF(data, filename, att_file=att_filename, verb)")
        return status

    keys = None
    nkeys = 0
    if isinstance(data, dict):
        keys = data.keys()
        nkeys = len(keys)
        if nkeys != 6:
            print(f"ERROR - Passed data dictionary does NOT contains all info")
            return status
    else:
        print(f"ERROR - Passed data is NOT a dictionary")
        return status

    att_filename = att_file
    if att_filename is None:
        print(f"WARNING - Attributes file not specified")
    elif not os.path.exists(att_filename):
        print(f"WARNING - Attributes file does NOT exists")

    # 
    # 	Do initial survey of variables and nested structures
    # 	to verify limitation on dimensions of arrays and nested structures
    # 	
    # 	LIMITATIONS:  4 dimensions on arrays and 4 nested structures
    # 
    # 	Use internal name structure for tracking any nested structures
    #
    dimsKeys = {}
    varDims = {}
    nameDims = {}
    for key in keys:
        shape = data[key].shape
        varDims[key] = shape
        for dim in shape:
            if dim not in dimsKeys:
                for dname in dims:
                    val = dims[dname]
                    if val == dim:
                        nameDims[dname] = dim
                        dimsKeys[dim] = dname
                        break


    now = datetime.datetime.now()
    compress = 'zlib'
    ncdf = None
    try:
        ncdf = netCDF4.Dataset(filename, 'w')#, format=netFormat)
    except:
        print(f"ERROR - Unable to open file {filename} for writing")
        return BAD_FILE


    ncdf.CreationTime = now.strftime("%Y-%m-%d %H:%M:%S.%f")
    uDims = {}
    for dimName in nameDims:
        uDims[dimName] = ncdf.createDimension(dimName, nameDims[dimName])

    for varName in keys:
        shape = list(varDims[varName])
        aDims = []
        for d in shape:
            aDims.append(uDims[dimsKeys[d]])
        ndim = len(shape)
        values = data[varName]
        atype = str(values.dtype)
        match atype:
            case 'float64':
                var = ncdf.createVariable(varName, (np.float64), (aDims))
            case 'float32':
                var = ncdf.createVariable(varName, (np.float32), (aDims))
            case 'int32':
                var = ncdf.createVariable(varName, (np.int32), (aDims))
            case 'int16':
                var = ncdf.createVariable(varName, (np.int16), (aDims))
            case 'int8':
                var = ncdf.createVariable(varName, (np.int8), (aDims))
        if ndim == 1:
            var[:] = values
        elif ndim == 2:
            var[:,:] = values
        elif ndim == 3:
            var[:,:,:] = values
        elif ndim == 4:
            var[:, :, :, :] = values
        else:
            print(f'ERROR - UNbale to handle variable {varName}, which has {ndim} dimensions')


    ncdf.close()

###############################################################################

def getFile(remurl, fout, key=None, usr=None):
    """
    Retrieves file associated with the remote URL link

    Parameters
    ----------
    iURL      Remote URL link for file to be retrieved
    fout      File name used to same remote file
    key         Site password, used if different than None
    usr         Site user name used if different than None


    Returns
    -------
    None
    """

    gcontext = None

    if re.search(r'https', remurl) is not None:
        gcontext = ssl.SSLContext(ssl.PROTOCOL_TLS)

    if usr is not None and key is not None:
        # create a password manager
        password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()

        # Add the username and password.
        # If we knew the realm, we could use it instead of None.
        password_mgr.add_password(None, remurl, usr, key)

        handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
        # create "opener" (OpenerDirector instance)
        opener = urllib.request.build_opener(handler)

        # Install the opener.
        # Now all calls to urllib.request.urlopen use our opener.
        urllib.request.install_opener(opener)

    try:
        with urllib.request.urlopen(remurl, context=gcontext) as response, open(fout, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
    except:
        msg = f"ERROR - Unable to get content from {remurl}\n"
        print(msg)
        raise


def retrieveSAVfile(archive, ydoy, remURL, verb):

    syear = str(ydoy.year)
    strydoy = ydoy.strftime("%Y%j")

    outdir = os.path.join(archive, syear)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    # URL = 'https://lasp.colorado.edu/eve/data_access/eve_data/fism/flare_hr_data/'+syear+'/'
    URL = remURL + syear + '/'

    remContent = (getContent(URL, verb=verb)).split("\n")

    pattern = re.compile('FISM_60sec_' + strydoy)
    files = []
    for line in remContent:
        if re.search(pattern, line):
            filename = (((line.split('href="'))[1]).split('<')[0]).split('>')[-1]
            files.append(filename)
            pass

    files.sort()
    fname = files[-1]
    remfile = URL + '/' + fname
    fout = os.path.join(outdir, fname)

    getFile(remfile, fout)

    return fout

def getContent(remurl, key=None, usr=None, verb=0):
    # print('Passed URL: '+URL+"\n"+str(type(URL)))

    """
    This function retrieves the content of a remote URL. It makes use of pycurl.

    Parameters
    ----------
    remurl      Remote URL
    key         Site password, used if different than None
    usr         Site user name used if different than None
    verb        Verbosity level

    Returns
    -------
    content     Remote content
    """
    content = None
    gcontext = None

    if re.search(r'https', remurl) is not None:
        gcontext = ssl.SSLContext(ssl.PROTOCOL_TLS)

    if usr is not None and key is not None:
        # create a password manager
        password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()

        # Add the username and password.
        # If we knew the realm, we could use it instead of None.
        password_mgr.add_password(None, remurl, usr, key)

        handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
        # create "opener" (OpenerDirector instance)
        opener = urllib.request.build_opener(handler)

        # Install the opener.
        # Now all calls to urllib.request.urlopen use our opener.
        urllib.request.install_opener(opener)

    try:
        f = urllib.request.urlopen(remurl, context=gcontext)
    except:

        try:
            print(urllib.error.HTTPError)
            print(urllib.error.HTTPError.reason())
        except:
            print("ERROR - Unable to retrieve content from " + remurl)
            return None

        return None

    content = f.read().decode('utf-8')

    return content


def writeFISMNetCDF(outfile, indict, verb):
    if os.path.exists(outfile):
        os.remove(outfile)

    if verb:
        print(f"INFO - Writing file {outfile}")


    now = datetime.datetime.now()
    compress = 'zlib'
    ncdf = netCDF4.Dataset(outfile, 'w')#, format=netFormat)

    # global attributes:
    ncdf.Year_DOY = indict['ydoy']
    ncdf.SecondOfDay = indict['sod']
    ncdf.CreationTime = now.strftime("%Y-%m-%d %H:%M:%S.%f")

    latDim = ncdf.createDimension('nLat', len(indict['lat']))
    lonDim = ncdf.createDimension('nLon', len(indict['lon']))


    if 'msk_irradiance' in indict:
        waveDim = ncdf.createDimension('nWave', len(indict['wavelength']))

        varlat = ncdf.createVariable('Lat', (np.float32), (latDim), zlib=1, complevel=5)
        varlon = ncdf.createVariable('Lon', (np.float32), (lonDim), zlib=1, complevel=5)
        varwav = ncdf.createVariable('Wavelengths', (np.float32), (waveDim), zlib=1, complevel=5)
        varirr = ncdf.createVariable('msk_irradiance', (np.float32), (latDim, lonDim, waveDim), zlib=1, complevel=5)
        varsza = ncdf.createVariable('sza', (np.float32), (latDim, lonDim), zlib=1, complevel=5)
        varise = ncdf.createVariable('iseclipse', (np.float32), (latDim, lonDim), zlib=1, complevel=5)
        varfis = ncdf.createVariable('fism2', (np.float32), (waveDim), zlib=1, complevel=5)

        varlat[:] = indict['lat']
        varlon[:] = indict['lon']
        varwav[:] = indict['wavelength']
        varirr[:,:,:] = indict['msk_irradiance']
        varsza[:,:] = indict['sza']
        varise[:,:] = indict['iseclipse']
        varfis[:] = indict['fism2']


    elif 'msk_irradiance_gitm' in indict:
        waveDim = ncdf.createDimension('nWave', len(indict['start_wv']))
        waveDim2 = ncdf.createDimension('nWave2', len(indict['fism2']))

        varlat = ncdf.createVariable('Lat', (np.float32), (latDim), zlib=1, complevel=5)
        varlon = ncdf.createVariable('Lon', (np.float32), (lonDim), zlib=1, complevel=5)
        varstartW = ncdf.createVariable('Start_wavelength', (np.float32), (waveDim), zlib=1, complevel=5)
        varendW = ncdf.createVariable('End_wavelength', (np.float32), (waveDim), zlib=1, complevel=5)
        varirr = ncdf.createVariable('msk_irradiance_gitm', (np.float32), (latDim, lonDim, waveDim), zlib=1, complevel=5)
        varsza = ncdf.createVariable('sza', (np.float32), (latDim, lonDim), zlib=1, complevel=5)
        varise = ncdf.createVariable('iseclipse', (np.float32), (latDim, lonDim), zlib=1, complevel=5)
        varfis = ncdf.createVariable('fism2', (np.float32), (waveDim2), zlib=1, complevel=5)
        varwave = ncdf.createVariable('Wavelengths', (np.float32), (waveDim2), zlib=1, complevel=5)


        varlat[:] = indict['lat']
        varlon[:] = indict['lon']
        varstartW[:] = indict['start_wv']
        varendW[:] = indict['end_wv']
        varirr[:,:,:] = indict['msk_irradiance_gitm']
        varsza[:,:] = indict['sza']
        varise[:,:] = indict['iseclipse']
        varfis[:] = indict['fism2']
        varwave[:] = indict['wavelength']


    else:
        print("ERROR - The passed dictionary with all data to be written to a file is unknown.")
        print(f"       Failed to write file {outfile}")

    ncdf.close()











