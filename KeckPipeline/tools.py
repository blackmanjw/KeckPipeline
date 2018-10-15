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


def swarp(files, output='output.fits', celestial_type='PIXEL'):
    """
        Run SWARP!

        Parameters
        ----------
        files: list
            A list of files to swarp.
        output: str
            Path and filename for the combined output *.fits file.
        celestial_type: str
            Define celetial type as in the swarp config/config.swarp file.
            Options are NAITVE, PIXEL, EQUATORIAL, GALACTIC, ECLIPTIC or SUPERGALACTIC.
    """
    kwargs = {
        'code': 'SWarp',
        'config': {
            'SUBTRACT_BACK': 'N',
            'IMAGEOUT_NAME': output,
            'RESCALE_WEIGHTS': 'N',
            'RESAMPLE': 'N',
            'CELESTIAL_TYPE': celestial_type,
            'INTERPOLATE': 'N',
            'BLANK_BADPIXELS': 'N',
            'COPY_KEYWORDS': 'OBJECT,CAMNAME,FWINAME,ITIME,DATE-OBS,FLSPECTR,HISTORY'
        },
        'temp_path': '.',
        'config_file': 'config/config.swarp'
    }

    swarp = aw.api.Astromatic(**kwargs)
    swarp.run(files)


def sex(files, output='output.cat', phot_aperture=10, back_size=256):
    """
        Run SExtractor!

        Parameters
        ----------
        files: list
            A list of files to SExtractor.
        output: str
            Path and filename for the catalog output  *.cat file.
        phot_aperture:
            Magnitude aperture diameter in pixels.
        back_size:
            Background mesh size in pixels.
    """
    kwargs = {
        'code': 'SExtractor',
        'config': {
            'CATALOG_NAME': output,
            'FILTER': 'N',
            'CLEAN': 'Y',
            'MASK_TYPE': 'CORRECT',
            'PHOT_APERTURES': str(phot_aperture),
            'DETECT_MINAREA': str(7),
            'DETECT_THRESH': str(7),
            'ANALYSIS_THRESH': str(7),
            'BACK_SIZE': str(back_size),
        },
        # Output parameters
        'params': ['NUMBER', 'XWIN_IMAGE', 'YWIN_IMAGE', 'MAG_AUTO', 'FLUX_APER', 'FLUXERR_APER'],
        'temp_path': '.',
        'config_file': 'config/config.sex'
    }

    sex = aw.api.Astromatic(**kwargs)
    sex.run(files)


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

def darkcombine(darks_dir='Darks/'):
    """
        This function combines Dark frames for each date and exposure time as per the filestructure created by "Rename".

        Parameters
        ----------
        darks_dir: str
            Define the folder where the dark files are located.

        Returns
        -------
        ????
    """

    # Make list of all Dark Directories
    darkdir = glob(darks_dir + '*/')
    print(darkdir)

    ## For each subdirectory in darkdir, combine the Dark Files (ATM only 30sec)
    for d in darkdir:
        keys = ['OBJECT', 'CAMNAME', 'FWINAME', 'ITIME', 'DATE-OBS']
        images = ImageFileCollection(d, keywords=keys)

        matches5 = (images.summary['ITIME'] < 6)
        dark5 = [d + x for x in images.summary['file'][matches5].tolist()]
        if dark5:
            print(dark5[2])
            swarp(dark5, output=darks_dir + 'Dark5sec' + d.split("-")[1] + d.split("-")[2][:-1] + '.fits')

        matches10 = (images.summary['ITIME'] == 10)
        dark10 = [d + x for x in images.summary['file'][matches10].tolist()]
        if dark10:
            swarp(dark10, output=darks_dir + 'Dark10sec' + d.split("-")[1] + d.split("-")[2][:-1] + '.fits')

        matches15 = (images.summary['ITIME'] == 15)
        dark15 = [d + x for x in images.summary['file'][matches15].tolist()]
        if dark15:
            swarp(dark15, output=darks_dir + 'Dark15sec' + d.split("-")[1] + d.split("-")[2][:-1] + '.fits')

        matches30 = (images.summary['ITIME'] == 30)
        dark30 = [d + x for x in images.summary['file'][matches30].tolist()]
        if dark30:
            swarp(dark30, output=darks_dir + 'Dark30sec' + d.split("-")[1] + d.split("-")[2][:-1] + '.fits')

        matches60 = (images.summary['ITIME'] == 60)
        dark60 = [d + x for x in images.summary['file'][matches60].tolist()]
        if dark60:
            swarp(dark60, output=darks_dir + 'Dark60sec' + d.split("-")[1] + d.split("-")[2][:-1] + '.fits')
        # print(y)


def darksubtract(dir='Flats/*', master_dark='Darks/Dark60sec0807.fits'):
    """
        This function subtracts the darks from files in a give directory.

        Parameters
        ----------
        ????

        Returns
        -------
        ????
    """
    # master_dark = input("Which Master Dark in the /Darks/ folder (or otherwise) would you like to use?")

    # if answer == 'Y':
    #    print("The script is now running....")
    # else:
    #    print("You have chosen to quit this program")
    #    raise SystemExit

    mdark = CCDData.read(master_dark, unit="adu")

    for d in glob(dir):
        keys = ['OBJECT', 'CAMNAME', 'FWINAME', 'ITIME', 'DATE-OBS', 'FLSPECTR', 'HISTORY']
        images = ImageFileCollection(d, keywords=keys, glob_exclude='d*', glob_include='*.fits')

        directory = d + '/dark_subtracted'
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Read in all files from /FLATS/ subdirectories and subtract the master_dark. The output is stored in 'dflat'.
        for flat, fname in images.hdus(return_fname=True):
            meta = flat.header
            meta['FILENAME'] = fname
            flat_exposure = flat.header['ITIME']
            flats = CCDData(data=flat.data.astype('float32'), meta=meta, unit="adu")
            dflat = (ccdproc.subtract_dark(flats, mdark, exposure_time='ITIME',
                                           exposure_unit=u.second,
                                           add_keyword={'HISTORY': 'Dark Subtracted'},
                                           scale=True))
            dflat.write(directory + '/d' + fname, overwrite=True)


def flatcombine(dir='Flats/*/dark_subtracted/'):
    """
        This function subtracts the darks from files in a give directory.

        Parameters
        ----------
        ????

        Returns
        -------
        ????
    """

    for d in glob(dir):

        directory = "/".join(d.split('/')[0:2]) + '/swarped'
        if not os.path.exists(directory):
            os.makedirs(directory)

        keys = ['OBJECT', 'CAMNAME', 'FWINAME', 'ITIME', 'DATE-OBS', 'FLSPECTR', 'HISTORY']
        images = ImageFileCollection(d, keywords=keys, glob_include='d*.fits')

        swarpfilter(d, dir, directory, images, keys, filter='H', lamp='on', camera='narrow', done='Dark Subtracted',
                    output='cKNarrowLampOnH')
        swarpfilter(d, dir, directory, images, keys, filter='H', lamp='off', camera='narrow', done='Dark Subtracted',
                    output='cKNarrowLampOffH')
        swarpfilter(d, dir, directory, images, keys, filter='H', lamp='on', camera='wide', done='Dark Subtracted',
                    output='cKWideLampOnH')
        swarpfilter(d, dir, directory, images, keys, filter='H', lamp='off', camera='wide', done='Dark Subtracted',
                    output='cKWideLampOffH')
        swarpfilter(d, dir, directory, images, keys, filter='Ks', lamp='on', camera='narrow', done='Dark Subtracted',
                    output='cKNarrowLampOnKs')
        swarpfilter(d, dir, directory, images, keys, filter='Ks', lamp='off', camera='narrow', done='Dark Subtracted',
                    output='cKNarrowLampOffKs')
        swarpfilter(d, dir, directory, images, keys, filter='Ks', lamp='on', camera='wide', done='Dark Subtracted',
                    output='cKWideLampOnKs')
        swarpfilter(d, dir, directory, images, keys, filter='Ks', lamp='off', camera='wide', done='Dark Subtracted',
                    output='cKWideLampOffKs')
        swarpfilter(d, dir, directory, images, keys, filter='J', lamp='on', camera='narrow', done='Dark Subtracted',
                    output='cNarrowLampOnJ')
        swarpfilter(d, dir, directory, images, keys, filter='J', lamp='off', camera='narrow', done='Dark Subtracted',
                    output='cKNarrowLampOffJ')
        swarpfilter(d, dir, directory, images, keys, filter='J', lamp='on', camera='wide', done='Dark Subtracted',
                    output='cKWideLampOnJ')
        swarpfilter(d, dir, directory, images, keys, filter='J', lamp='off', camera='wide', done='Dark Subtracted',
                    output='cKWideLampOffJ')


# ------------------------------ REQUIRED FUNCTIONS -----------------------------

def swarpfilter(d, dir, directory, images, keys, filter, lamp, camera, done, output):
    """
        This function runs a swarp on a subset of images from an ImageFileCollection according to the chosen parameters
        (Filter, Band, Camera etc.) This function is only used in other functions like flatcombine.
    """
    filt = images.files_filtered(FWINAME=filter, FLSPECTR=lamp, CAMNAME=camera, HISTORY=done)
    files = [d + x for x in filt.tolist()]
    print(files)
    if files:
        swarp(files, output=directory + '/' + output + '.fits')
