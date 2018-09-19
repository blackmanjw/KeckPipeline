import os
from astropy.io import fits
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

## Rename Files, Move and Create Folders

def movedata(startdate, numdays):  # where startdate is in the form '01-01-2018' and numdays is just an integer.

    # Creating a list of dates starting from the startdate, and lasts a number of days that you input.
    datelist = pd.date_range(startdate, periods=numdays, freq='D').strftime('%d-%m-%Y').tolist()

    # Create folders named after the list of dates.
    for folder in datelist:
        if not os.path.exists('datelist'):
            os.mkdir(os.path.join('Data', str(folder)))


def rename(dir):
    data = []
    for file in os.listdir("./" + dir):  # put in your path directory
        if file.endswith(".fits"):  # what does the file end with?
            data.append(os.path.join(dir, file))

    n = len(data)

    obj = []
    itime = []
    filt = []
    renamed = []
    Date_1 = []
    count = []
    for i in range(0, n):
        header = fits.getheader(data[i])
        obj.append(header['OBJECT'])
        itime.append(header['ITIME'])
        filt.append(header['FWINAME'])
        Name, Date, Number, Ext = data[i].split(".")
        Date_1.append(datetime.strptime(Date, "%Y%m%d").date())
        if obj[i] in obj:
            count = obj.count(obj[i])
        if ('Lamp' in obj[i] or 'Flat' in obj[i]):
            renamed.append(('Flats/' + str(Date_1[i]) + '/' + 'K' + header['OBJECT'] + '_' + str(count) + ".fits"))
            os.makedirs(os.path.dirname('Flats/' + str(Date_1[i]) + '/'), exist_ok=True)
        elif 'dark' in obj[i]:
            renamed.append(('Darks/' + str(Date_1[i]) + '/' + 'K' + header['OBJECT'] + header['FWINAME'] + '_' + str(
                count) + ".fits"))
            os.makedirs(os.path.dirname('Darks' + str(Date_1[i]) + '/'), exist_ok=True)
        elif 'sky' in obj[i]:
            renamed.append(('Skys/' + str(Date_1[i]) + '/' + 'K' + header['OBJECT'] + header['FWINAME'] + '_' + str(
                count) + ".fits"))
            os.makedirs(os.path.dirname('Skys/' + str(Date_1[i]) + '/'), exist_ok=True)
        else:
            renamed.append(('Data/' + str(Date_1[i]) + '/' + 'K' + header['OBJECT'] + header['FWINAME'] + '_' + str(
                count) + ".fits"))
            os.makedirs(os.path.dirname('Data/' + str(Date_1[i]) + '/'), exist_ok=True)
        os.rename(data[i], renamed[i])

    # Name,Date,Number,Ext=line.split(".")
    lists = [data, obj, Date_1, itime, filt, renamed]
    data_headers = pd.concat([pd.Series(x) for x in lists], axis=1)
    print(data_headers)

    return data_headers, obj


## Create *.list file of darkframes

def darklist(dir):

    darkfiles = []
    for file in os.listdir("./" + dir):
        if file.endswith(".fits"):
            darkfiles.append(os.path.join(dir, file))

    with open(dir + '/dark.list', 'w') as f:
        for item in darkfiles:
            f.write("%s\n" % item)
    return darkfiles

### Combine Dark Frames

def darkcombine():

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