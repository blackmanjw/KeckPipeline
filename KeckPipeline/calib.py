import os
from astropy.io import fits
import pandas as pd
import numpy as np
from datetime import datetime, timedelta


def movedata(startdate, numdays):
    data = []
    for file in os.listdir("./2018"):
        if file.endswith(".fits"):
            data.append(os.path.join("2018", file))

    datelist = pd.date_range(startdate, periods=numdays, freq='D').strftime('%d-%m-%Y').tolist()

    for folder in datelist:
        os.mkdir(os.path.join('2018', str(folder)))

    n = len(data)

    obj = []
    itime = []
    filt = []
    renamed = []
    Date_1 = []
    for i in range(0, n):
        header = fits.getheader(data[i])
        obj.append(header['OBJECT'])
        itime.append(header['ITIME'])
        filt.append(header['FWINAME'])
        Name, Date, Number, Ext = data[i].split(".")
        Date_1.append(datetime.strptime(Date, "%Y%m%d").date())
        renamed.append(('2018/' + 'K' + header['OBJECT'] + header['FWINAME'] + '' + ".fits"))

    # Name,Date,Number,Ext=line.split(".")
    lists = [data, obj, Date_1, itime, filt, renamed]
    data_headers = pd.concat([pd.Series(x) for x in lists], axis=1)

    return Date, data_headers


## CREATE MASTER DARK FRAME

def darks():
    ## Create List of Darks

    darkfiles = []
    for file in os.listdir("./Darks"):
        if file.endswith(".fits"):
            darkfiles.append(os.path.join("Darks", file))

    with open('Darks/dark.list', 'w') as f:
        for item in darkfiles:
            f.write("%s\n" % item)

    # Import list of Dark Frames

    darklist = [line.rstrip('\n') for line in open('Darks/dark.list')]

    print("\n" 'DARK Source Files:')

    print("\n".join(darkfiles))

    print("\n" 'DARK Exposure Times:')

    darksexp = []
    for line in darklist:
        header = fits.getheader(line)
        darksexp.append(header['ITIME'])
        print(header['ITIME'])

    darks = []
    for file in darklist:
        darks.append(fits.getdata(file))

    mediandark = np.median(darks, axis=0)

    header = fits.getheader(darklist[0])

    header['HISTORY'] = 'Median combined'

    fits.writeto('dark_300718_30s.fits', mediandark, header, overwrite=True)

    return mediandark


## SUBTRACT DARK FRAMES FROM DOME FLATS

## Create List of Dome Flats

def flats():
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