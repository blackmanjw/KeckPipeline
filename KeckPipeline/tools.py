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

import pandas as pd
import numpy as np
import shutil
from datetime import datetime, timedelta
from glob import glob
from ccdproc import ImageFileCollection
import warnings
from astropy.io import fits
import astropy.units as u
from astropy.nddata import CCDData
import ccdproc
from itertools import product
import os
from astropy.io.fits import getdata
import astromatic_wrapper as aw
import numpy as np
from distutils import dir_util
from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats

warnings.filterwarnings("ignore")


# ------------------------------ USEFUL FUNCTIONS -----------------------------

def keep_going(text="Do you wish to continue? Answer Y or N."):
    """
        This function asks the user whether or not they want the script to proceed.

        The user needs to type 'Y' (no quotes), else the script will abort. Hmm.

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


def swarp(files, output='output.fits', celestial_type='EQUATORIAL'):
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
            'COMBINE_TYPE': 'MEDIAN',
            'PIXEL_SCALE': '0',
            'CELESTIAL_TYPE': celestial_type,
            'INTERPOLATE': 'N',
            'BLANK_BADPIXELS': 'N',
            'COPY_KEYWORDS': 'OBJECT,CAMNAME,FWINAME,ITIME,DATE-OBS,OBSDATE,FLSPECTR,HISTORY'
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
        'params': ['NUMBER', 'BACKGROUND', 'XWIN_IMAGE', 'YWIN_IMAGE', 'MAG_AUTO', 'FLUX_APER', 'FLUXERR_APER'],
        'temp_path': '.',
        'config_file': 'config/config.sex'
    }

    sex = aw.api.Astromatic(**kwargs)
    sex.run(files)


from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats


def backgroundestimate(file, size=10):
    """
        This function estimates the background of a fits images via iteration and masking of any sources.

        Required imports: from photutils import make_source_mask
                          from astropy.stats import sigma_clipped_stats

        Parameters
        ----------
        file:
            Path to fits file.
        output:
            The mean background count.
    """
    data = CCDData.read(file, unit="adu")
    mask = make_source_mask(data, snr=2, npixels=5, dilate_size=size)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    print(mean)
    return mean


# ------------------------------ PIPELINE FUNCTIONS -----------------------------

def rename(source_dir, dest_dir):
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
    keep_going(
        text="This script will backup the original folder to dest_dir/Source/** and remove the original folder. It will make copies of the  original files and rename them in directories called Darks, Flats, etc. Do you wish to continue? Answer Y or N.")

    ## Backup Original Source Folder
    dir_util.copy_tree(source_dir, dest_dir + '/Source')

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
            renamed.append(
                (dest_dir + '/Flats/' + str(datemod[i]) + '/' + 'K' + header['OBJECT'] + str(count) + ".fits"))
            os.makedirs(os.path.dirname(dest_dir + '/Flats/' + str(datemod[i]) + '/'), exist_ok=True)
        elif ('Dark' in obj[i]) or ('dark' in obj[i]):
            renamed.append(
                (dest_dir + '/Darks/' + str(datemod[i]) + '/' + 'K' + header['OBJECT'] + str(count) + ".fits"))
            os.makedirs(os.path.dirname(dest_dir + '/Darks/' + str(datemod[i]) + '/'), exist_ok=True)
        elif ('Sky' in obj[i]) or ('sky' in obj[i]):
            renamed.append(
                (dest_dir + '/Skys/' + str(datemod[i]) + '/' + 'K' + header['OBJECT'] + header['FWINAME'] + str(
                    count) + ".fits"))
            os.makedirs(os.path.dirname(dest_dir + '/Skys/' + str(datemod[i]) + '/'), exist_ok=True)
        else:
            renamed.append((dest_dir + '/Objects/' + header['OBJECT'].upper() + '/' + str(datemod[i]) + '/' + 'K' +
                            header['OBJECT'].upper() +
                            header['FWINAME'] + str(
                        count) + ".fits"))
            os.makedirs(
                os.path.dirname(dest_dir + '/Objects/' + header['OBJECT'].upper() + '/' + str(datemod[i]) + '/'),
                exist_ok=True)
        os.rename(data[i], renamed[i])

    ## REMOVE LEFT OVER original Folders
    shutil.rmtree(source_dir)

    # Name,Date,Number,Ext=line.split(".")
    lists = [data, mod, datemod, itime, flatmod, renamed]
    data_headers = pd.concat([pd.Series(x) for x in lists], axis=1)
    # print(data_headers)

    return data_headers


# -----------------------------

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


# -----------------------------

def darksubtract(dir='Flats', master_dark='Darks/Dark60sec0807.fits'):
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

    if dir == "Flats":
        dir = 'Flats/*/'
    elif dir == "Skys":
        dir = 'Skys/*/'
    elif dir == "Objects":
        dir = 'Objects/*/*/'

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
                                           add_keyword={'HISTORY': 'Dark Subtracted',
                                                        'OBSDATE': flat.header['DATE-OBS']},
                                           scale=True))
            dflat.write(directory + '/d' + fname, overwrite=True)


# -----------------------------

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

        keys = ['OBJECT', 'CAMNAME', 'FWINAME', 'ITIME', 'OBSDATE', 'FLSPECTR', 'HISTORY']
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


# -----------------------------

def flatprocess(dir='Flats/*/swarped/'):
    """

        Parameters
        ----------
        ????

        Returns
        -------
        ????
    """
    print(glob(dir))
    for d in glob(dir):

        keys = ['OBJECT', 'CAMNAME', 'FWINAME', 'ITIME', 'OBSDATE', 'FLSPECTR', 'HISTORY']
        images = ImageFileCollection(d, keywords=keys, glob_include='c*.fits')
        cameras = ['wide', 'narrow']
        filters = ['Ks', 'H', 'J']
        flats = []

        for camera, filter in product(cameras, filters):

            for hdu, fname in images.hdus(CAMNAME=camera, FWINAME=filter, return_fname=True):
                meta = hdu.header
                meta['filename'] = fname
                flats.append(ccdproc.CCDData(data=hdu.data.astype('float32'), meta=meta, unit="adu"))
            if flats:
                MasterFlat = (ccdproc.subtract_dark(flats[1], flats[0], exposure_time='ITIME',
                                                    exposure_unit=u.second,
                                                    add_keyword={'HISTORY2': 'Lamp Off Subtracted'},
                                                    scale=False))

                # Normalize Flat
                MasterFlat.data = MasterFlat / np.ma.average(MasterFlat)
                MasterFlat.write(
                    'Flats/KFlat' + MasterFlat.meta['CAMNAME'].title() + MasterFlat.meta['FWINAME'].title() +
                    MasterFlat.meta['OBSDATE'].split("-")[1] + MasterFlat.meta['OBSDATE'].split("-")[2] + '.fits',
                    overwrite=True)
                flats[:] = []


# -----------------------------

def flatcorrect(dir='Skys', flatdate='2018-08-07'):
    """

        Parameters
        ----------
        ????

        Returns
        -------
        ????
    """

    if dir == "Skys":
        dir = 'Skys/*/dark_subtracted/'
    elif dir == "Objects":
        dir = 'Objects/*/*/dark_subtracted/'

    keys = ['OBJECT', 'CAMNAME', 'FWINAME', 'ITIME', 'OBSDATE', 'FLSPECTR', 'HISTORY']
    flatfiles = ImageFileCollection('Flats', keywords=keys)
    cameras = ['wide', 'narrow']
    filters = ['Ks', 'H', 'J']
    add_filters = {'OBSDATE': flatdate}
    flats = {}

    ## IMPORT MASTER FLATS for chosen date in dictionary flats={}

    for camera, filter in product(cameras, filters):

        for hdu, fname in flatfiles.hdus(CAMNAME=camera, FWINAME=filter, return_fname=True, **add_filters):
            meta = hdu.header
            meta['filename'] = fname
            flats[fname[:-5]] = ccdproc.CCDData(data=hdu.data.astype('float32'), meta=meta, unit="adu")

    ## APPLY FLAT CORRECTION TO IMAGES

    for d in glob(dir):

        if dir == "Skys/*/dark_subtracted/":
            directory = "/".join(d.split('/')[0:2]) + '/flat_corrected'
        elif dir == "Objects/*/*/dark_subtracted/":
            directory = "/".join(d.split('/')[0:3]) + '/flat_corrected'

        if not os.path.exists(directory):
            os.makedirs(directory)

        images = ImageFileCollection(d, keywords=keys, glob_include='d*.fits')
        skys = []

        for camera, filter in product(cameras, filters):

            for hdu, fname in images.hdus(CAMNAME=camera, FWINAME=filter, return_fname=True):
                meta = hdu.header
                meta['filename'] = fname
                sky = ccdproc.CCDData(data=hdu.data.astype('float32'), meta=meta, unit="adu")

                if (meta['CAMNAME'] == 'wide') and (meta['FWINAME'] == 'Ks'):
                    cflat = ccdproc.flat_correct(sky,
                                                 flats['KFlatWideKs' + flatdate.split('-')[1] + flatdate.split('-')[2]])
                    cflat.write(directory + '/f' + fname, overwrite=True)
                elif (meta['CAMNAME'] == 'narrow') and (meta['FWINAME'] == 'Ks'):
                    cflat = ccdproc.flat_correct(sky, flats[
                        'KFlatNarrowKs' + flatdate.split('-')[1] + flatdate.split('-')[2]])
                    cflat.write(directory + '/f' + fname, overwrite=True)
                elif (meta['CAMNAME'] == 'wide') and (meta['FWINAME'] == 'H'):
                    cflat = ccdproc.flat_correct(sky,
                                                 flats['KFlatWideH' + flatdate.split('-')[1] + flatdate.split('-')[2]])
                    cflat.write(directory + '/f' + fname, overwrite=True)
                elif (meta['CAMNAME'] == 'narrow') and (meta['FWINAME'] == 'H'):
                    cflat = ccdproc.flat_correct(sky, flats[
                        'KFlatNarrowH' + flatdate.split('-')[1] + flatdate.split('-')[2]])
                    cflat.write(directory + '/f' + fname, overwrite=True)
                elif (meta['CAMNAME'] == 'wide') and (meta['FWINAME'] == 'J'):
                    cflat = ccdproc.flat_correct(sky,
                                                 flats['KFlatWideJ' + flatdate.split('-')[1] + flatdate.split('-')[2]])
                    cflat.write(directory + '/f' + fname, overwrite=True)
                elif (meta['CAMNAME'] == 'narrow') and (meta['FWINAME'] == 'J'):
                    cflat = ccdproc.flat_correct(sky, flats[
                        'KFlatNarrowJ' + flatdate.split('-')[1] + flatdate.split('-')[2]])
                    cflat.write(directory + '/f' + fname, overwrite=True)


# -----------------------------

def skycombine(dir='Objects'):
    """
        This function combines the skys from files in a give directory, after they have been dark subtracted and flat
        corrected.

        Parameters
        ----------
        ????

        Returns
        -------
        ????
    """

    if dir == "Objects":
        dir = 'Objects/*/*/flat_corrected/'

    for d in glob(dir):

        directory = "/".join(d.split('/')[0:2]) + '/swarped'
        if not os.path.exists(directory):
            os.makedirs(directory)

        keys = ['OBJECTS', 'ITIME', 'FWINAME', 'OBSDATE', 'CAMNAME', 'HISTORY', 'FLSPECTR']
        images = ImageFileCollection(d, keywords=keys, glob_include='f*.fits')

        swarpfilter(d, dir, directory, images, keys, filter='H', lamp='*', camera='narrow', done='Dark Subtracted',
                    output='cKSkyNarrowH')
        swarpfilter(d, dir, directory, images, keys, filter='H', lamp='*', camera='wide', done='Dark Subtracted',
                    output='cKSkyWideH')
        swarpfilter(d, dir, directory, images, keys, filter='J', lamp='*', camera='narrow', done='Dark Subtracted',
                    output='cKSkyNarrowJ')
        swarpfilter(d, dir, directory, images, keys, filter='J', lamp='*', camera='wide', done='Dark Subtracted',
                    output='cKSkyWideJ')
        swarpfilter(d, dir, directory, images, keys, filter='Ks', lamp='*', camera='narrow', done='Dark Subtracted',
                    output='cKSkyNarrowKs')
        swarpfilter(d, dir, directory, images, keys, filter='Ks', lamp='*', camera='wide', done='Dark Subtracted',
                    output='cKSkyWideKs')
        swarpfilter(d, dir, directory, images, keys, filter='Lp', lamp='*', camera='narrow', done='Dark Subtracted',
                    output='cKSkyNarrowLp')
        swarpfilter(d, dir, directory, images, keys, filter='Lp', lamp='*', camera='wide', done='Dark Subtracted',
                    output='cKSkyNarrowLp')


# ---------------------------

def skycorrect(dir='Objects'):
    """
        This function applys the sky correction to the object files. This takes the form:

        OBJECT - SKY (mean object background / mean sky background)

        Parameters
        ----------
        dir: str
            Target directory containing subdirectories for each object.

        Returns
        -------
        ????
    """

    if dir == "Objects":
        dir = 'Objects/*/*/flat_corrected/'

    keys = ['OBJECT', 'CAMNAME', 'FWINAME', 'ITIME', 'OBSDATE', 'FLSPECTR', 'HISTORY']
    cameras = ['wide', 'narrow']
    filters = ['Ks', 'H', 'J']
    skys = {}
    objects = {}

    ## Import MASTER SKYS from swarped directory for all date folders under Skys/*

    for d in glob('Skys/*/swarped'):

        files = ImageFileCollection(d, keywords=keys, glob_include='*.fits')

        for camera, filter in product(cameras, filters):

            for hdu, fname in files.hdus(CAMNAME=camera, FWINAME=filter, return_fname=True):
                meta = hdu.header
                meta['filename'] = fname
                skys[fname[:-5] + meta['OBSDATE'].split("-")[1] + meta['OBSDATE'].split("-")[2]] = [
                    ccdproc.CCDData(data=hdu.data.astype('float32'), meta=meta, unit="adu"),
                    backgroundestimate(d + '/' + fname)]
                print(skys)
    ## APPLY SKY CORRECTION to Images


''' 
    for d in glob(dir):

        directory = "/".join(d.split('/')[0:3]) + '/sky_corrected'
        if not os.path.exists(directory):
            os.makedirs(directory)

        images = ImageFileCollection(d, keywords=keys, glob_include='*.fits')

        for camera, filter in product(cameras, filters):

            for hdu, fname in images.hdus(CAMNAME=camera, FWINAME=filter, return_fname=True):
                meta = hdu.header
                meta['filename'] = fname
                objects[fname[:-5]+ meta['OBSDATE'].split("-")[1] + meta['OBSDATE'].split("-")[2]] = [ccdproc.CCDData(data=hdu.data.astype('float32'), meta=meta, unit="adu"),backgroundestimate(d + '/' + fname)]
                print(objects)
                #MasterFlat.data = MasterFlat / np.ma.average(MasterFlat)
                #MasterFlat.write(
                #    'Flats/KFlat' + MasterFlat.meta['CAMNAME'].title() + MasterFlat.meta['FWINAME'].title() +
                 #   MasterFlat.meta['OBSDATE'].split("-")[1] + MasterFlat.meta['OBSDATE'].split("-")[2] + '.fits',
                 #   overwrite=True)
               # flats[:] = []

'''


# Normalize Flat
# MasterFlat.data = MasterFlat / np.ma.average(MasterFlat)

# backgroundestimate(file='Objects/OB171434/2018-08-07/KOB171434Ks1.fits')

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