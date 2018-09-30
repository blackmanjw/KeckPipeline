############################################################################
#                       KECK NIRC2 AO Reduction Tools                      #
#                                                                          #
#                             REVISION: 0.0.1                              #
#                                                                          #
#                             Joshua Blackman                              #
#                           Aikaterini Vandorou                            #
#                          Univeristy of Tasmania                          #
#                                                                          #
#                          Last Change: Oct 2018                           #
############################################################################
#
#
# This is a pipeline to reduce AO images from the NIRC2 imager on KECK II.
#
#

import os
from astropy.io import fits
import astromatic_wrapper as aw
import glob
import pandas as pd
import numpy as np
import shutil
from datetime import datetime, timedelta


# ------------------------------ USEFUL FUNCTIONS -----------------------------

def keep_going(text="Do you wish to continue? Answer Y or N."):
    """
        This function asks the user whether or not they want the script to proceed.

        The user needs to type 'Y' (no quotes), else the script will abort.

        Parameters
        ----------
        text: str
        Modify text to present to the user. Default is "Do you wish to continue? Answer Y or N."
        """
    answer = input(text)

    if answer == 'Y':
        print("The script is now running....")
    else:
        print("You have chosen to quit this program")
        raise SystemExit


# ------------------------------ PIPELINE FUNCTIONS -----------------------------

def rename(source_dir,dest_dir):
    """
        This function takes the downloaded files from the KOA website, renames them and sorts into folders named Objects, Darks, Flats and Skys.
        Objects are sorted by object name (read from the fits header), then date. Calibration files are sorted by date.
        This function will DELETE the original folder but backup the folder into a directory called Source.
        Before proceeding the function will ask the user if they wish to proceed.

        As with the other functions in tools.py, this function is run via calib.py.

        This function can be run with "python3 calib.py -r -s <source_dir> <dest_dir>" where <soure_folder> is the path to your data you wish to move and rename.

        Parameters
        ----------
        source_dir: str
            Define the folder where the source files are located.
        dest_dir: str
            Define the output folder. The Object, Darks etc. folders will be created under this parent.

        Returns
        -------
        data_headers
            Table with columns showing input files, output files, object name and exposure time.
    """
    keep_going(text="This script will backup the original folder to dest_dir/Source/** and remove the original folder. It will make copies of the  original files and rename them in directories called Darks, Flats, etc. Do you wish to continue? Answer Y or N.")

    ## Backup Original Source Folder
    shutil.copytree(source_dir, dest_dir + '/Source')

    data = []
    for file in os.listdir("./" + source_dir):  # put in your path directory
        if file.endswith(".fits"):  # what does the file end with?
            data.append(os.path.join(source_dir, file))

    n = len(data)
    obj, itime, filt, renamed, datemod, count, flatmod, mod = ([] for i in range(8))
    for i in range(0, n):
        header = fits.getheader(data[i])
        Name, Date, Number, Ext = data[i].split(".")
        obj.append(header['OBJECT'])
        itime.append(header['ITIME'])
        filt.append(header['FWINAME'])
        mod.append((header['OBJECT'] + header['FWINAME']))
        flatmod.append((header['OBJECT'] + header['FWINAME'] + Date))
        datemod.append(datetime.strptime(Date, "%Y%m%d").date())
        if flatmod[i] in flatmod:
            count = flatmod.count(flatmod[i])
        if ('Lamp' in obj[i] or 'Flat' in obj[i]):
            renamed.append((dest_dir + '/Flats/' + str(datemod[i]) + '/' + 'K' + header['OBJECT'] + str(count) + ".fits"))
            os.makedirs(os.path.dirname(dest_dir + '/Flats/' + str(datemod[i]) + '/'), exist_ok=True)
        elif ('Dark' in obj[i]) or ('dark' in obj[i]):
            renamed.append((dest_dir + '/Darks/' + str(datemod[i]) + '/' + 'K' + header['OBJECT'] + str(count) + ".fits"))
            os.makedirs(os.path.dirname(dest_dir + '/Darks/' + str(datemod[i]) + '/'), exist_ok=True)
        elif ('Sky' in obj[i]) or ('sky' in obj[i]):
            renamed.append((dest_dir + '/Skys/' + str(datemod[i]) + '/' + 'K' + header['OBJECT'] + header['FWINAME'] + str(
                count) + ".fits"))
            os.makedirs(os.path.dirname(dest_dir + '/Skys/' + str(datemod[i]) + '/'), exist_ok=True)
        else:
            renamed.append((dest_dir + '/Objects/' + header['OBJECT'].upper() + '/' + str(datemod[i]) + '/' + 'K' + header['OBJECT'].upper() +
                            header['FWINAME'] + str(
                        count) + ".fits"))
            os.makedirs(os.path.dirname(dest_dir + '/Objects/' + header['OBJECT'].upper() + '/' + str(datemod[i]) + '/'), exist_ok=True)
        os.rename(data[i], renamed[i])

    ## REMOVE LEFT OVER original Folders
    shutil.rmtree(source_dir)

    # Name,Date,Number,Ext=line.split(".")
    lists = [data, mod, datemod, itime, flatmod, renamed]
    data_headers = pd.concat([pd.Series(x) for x in lists], axis=1)
    print(data_headers)

    return data_headers

### Combine Dark Frames

def darkcombine():
    """
        LLALALALALALALALALALALALLALALLALA
    """
    darkfiles=darklist(dir)

    dark_list = [line.rstrip('\n') for line in open('Darks/dark.list')]

    print("\n" 'DARK Source Files:')

    print("\n".join(darkfiles))

    print("\n" 'DARK Exposure Times:')

    darksexp = []
    for line in dark_list:
        header = fits.getheader(line)
        darksexp.append(header['ITIME'])
        print(header['ITIME'])

    darks = []
    for file in dark_list:
        darks.append(fits.getdata(file))

    mediandark = np.median(darks, axis=0)

    header = fits.getheader(dark_list[0])

    header['HISTORY'] = 'Median combined'

    fits.writeto('dark_300718_30s.fits', mediandark, header, overwrite=True)

    return mediandark




## Create *.list file of darkframes

def swarp():

    data = []
    for file in os.listdir("./" + Darks):  # put in your path directory
        if file.endswith(".fits"):  # what does the file end with?
            data.append(os.path.join(Darks, file))
    print(data)
    kwargs = {
        'code': 'SWarp',
        'config': {
            'SUBTRACT_BACK': 'N',
            'IMAGEOUT_NAME': 'darks.fits',
            'RESCALE_WEIGHTS': 'N',
            'RESAMPLE': 'N',
            'INTERPOLATE': 'Y',
            'BLANK_BADPIXELS': 'N',
        },
        'temp_path': '.',
        'config_file': '../../config/config.swarp'
    }

    darks=glob.glob('../2018-08-07/K*30*.fits')
    print(darks)
    swarp = aw.api.Astromatic(**kwargs)
    #sextractor = aw.api.Astromatic(**kwargs)
    swarp.run(darks)

def darklist(dir):

    dir=[x[0] for x in os.walk('Darks')]

    n=len(dir)
    darkfiles = []
    for i in range (1,n):
        for file in os.listdir(dir[i]):
            if file.endswith(".fits"):
                darkfiles.append(os.path.join(dir[i], file))
        with open(dir[i] + '/dark.list', 'w') as f:
            for item in darkfiles:
                f.write("%s\n" % item)

    return darkfiles


## SUBTRACT DARK FRAMES FROM DOME FLATS

## Create List of Dome Flats

## Create *.list file of Dome Flats

def flatlist(dir):

    flats = []
    for file in os.listdir("./" + dir):
        if file.endswith(".fits"):
            flats.append(os.path.join(dir, file))

    with open(dir + '/flats.list', 'w') as f:
        for item in flats:
            f.write("%s\n" % item)
    return flats

def flatcombine():

    flats = flatlist()
    mediandark = darks()
    domeflatfiles = []
    for file in os.listdir("./DomeFlats"):
        if file.endswith(".fits"):
            domeflatfiles.append(os.path.join("DomeFlats", file))

    print("\n" 'DOME FLAT LampOFF Source Files:')
    print("\n".join(domeflatfiles))

    with open('DomeFlats/domeflat.list', 'w') as f:
        for item in domeflatfiles:
            f.write("%s\n" % item)

    print("\n" 'DOME FLAT Exposure Times:')

    for line in domeflatfiles:
        header = fits.getheader(line)
        print(header['ITIME'])

    if not os.path.exists('./Calib'):
        os.makedirs('./Calib')

    domelist = [line.rstrip('\n') for line in open('DomeFlats/domeflat.list')]

    domeflat_dark = []
    n = 0
    for line in domelist:
        Name, Date, Number, Ext = line.split(".")
        n = n + 1
        domeflat_dark.append('Calib/' + Name[10:-3] + '_' + Date + '_' + 'dark' + '_' + str(n) + ".fits")

    print("\n" 'DOME FLATS LampOFF (Bias Subtracted) Calib Output:')
    print("\n".join(domeflat_dark))

    with open('Calib/domeflat_dark.list', 'w') as f:
        for item in domeflat_dark:
            f.write("%s\n" % item)

    n = len(domeflatfiles)

    for i in range(0, n):
        data, header = fits.getdata(domeflatfiles[i], header=True)
        domeout = data - mediandark
        header['HISTORY'] = 'Bias subtracted'
        fits.writeto(domeflat_dark[i], domeout, header, overwrite=True)

    return domeflat_dark, domeflatfiles


## CREATE MASTER DOME FLATS

##Ks BAND

def domeflat():
    domeflat_dark, domeflatfiles = flats()
    Date, data_headers = movedata(startdate, numdays)
    ksflat_onstack = []
    ksflat_offstack = []

    for file in domeflat_dark:
        if "LampOn" in files:
            data, header = fits.getdata(file, header=True)
            ksflat_onstack.append(data)
        if "LampOff" in file:
            data, header = fits.getdata(file, header=True)
            ksflat_offstack.append(data)

    ksflat_on = np.median(ksflat_onstack, axis=0)
    ksflat_off = np.median(ksflat_offstack, axis=0)

    ksflat = ksflat_on - ksflat_off
    m = np.mean(ksflat)
    ksflat = ksflat / m

    print(ksflat)

    header['HISTORY'] = 'Combined and Normalized Flat Field'
    fits.writeto('Calib/ksflat' + '_' + Date + ".fits", ksflat, header, overwrite=True)

    return ksflat


## COMBINE SKY FLATS

