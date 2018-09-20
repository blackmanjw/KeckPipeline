# KECKPipeline

Authors : Joshua Blackman, joshua.blackman@utas.edu.au;
	      Aikaterini Vandorou, aikaterini.vandorou@utas.edu.au;
	      Jean-Phillipe Beaulieu,=

This is a pipeline for performing calibration corrections on NIRC2 images obtained on the Keck II telescope.

# Documentation and Installation

**Rename Files** 
```python
python3 calib.py -r --rename
```
**Combine Dark Frames**

# Naming Conventions

Running

```python
python3 calib.py -r 
```

with the -r (--rename) flag moves all raw data obtained from the KOA website (https://koa.ipac.caltech.edu/cgi-bin/KOA/nph-KOAlogin into a more convenient directory structure. That structure is shown below:

**Folder Structure**
    
    ├── Darks                  
    │   ├── 2018-08-05                
    │       ├── KDark30sec1.fits   
    │       └── KDark30sec2.fits  
    ├── Flats               
    │   ├── 2018-08-05                
    │       ├── KNarrowLampOffKs1.fits  
    │       ├── KNarrowLampOnKs1.fits
    │       ├── KNarrowLampOffH1.fits  
    │       ├── KNarrowLampOffH2.fits  
    │       ├── KNarrowLampOnH1.fits 
    │       ├── KNarrowLampOnH2.fits      
    │       ├── KWideLampOnKs1.fits
    │       └── KWideLampOnKs2.fits 
    ├── Objects               
    │   └── ob161455          
    │       ├── 2018-05-06 
    │           ├── KNarrowLampOnH1.fits 
    │           ├── KNarrowLampOnH2.fits 
    │           └── KWideLampOnKs2.fits          
    │       └── 2018-05-07  
    │           ├── KNarrowLampOnH1.fits 
    │           ├── KNarrowLampOnH2.fits 
    │           └── KWideLampOnKs2.fits 
    │   └── mb101255          
    │       └── 2018-05-06    
    │           ├── KNarrowLampOnH1.fits 
    │           └── KWideLampOnKs2.fits
    ├── Skys                  # Automated tests (alternatively `spec` or `tests`)
    │   ├── 2018-05-05        # Load and stress tests
    │   └── 2018-05-06        # Unit tests
    ├── Source                # Backup of Original KOA Download
    ├── calib.py
    └── tools.py
