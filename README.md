# KECKPipeline

Authors : Joshua Blackman, joshua.blackman@utas.edu.au;
	      Aikaterini Vandorou, aikaterini.vandorou@utas.edu.au;
	      Jean-Phillipe Beaulieu

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
python3 calib.py -r -s folder 
```

with the -r (--rename) flag moves all raw data obtained from the KOA website (https://koa.ipac.caltech.edu/cgi-bin/KOA/nph-KOAlogin) from a source folder into a more convenient directory structure. The source folder can be specificed using the '-s' flag. The resulting structure is shown below:

**Raw Data Folder Structure**
    
    ├── Darks 
    │   ├── 2018-08-05                
    │       ├── KDark30sec1.fits             # K Dark30sec 1 (Keck, Dark Type, Image Number)
    │       └── KDark30sec2.fits 
    ├── Flats               
    │   ├── 2018-08-05                
    │       ├── KNarrowLampOffKs1.fits       # K NarrowLampOff Ks 1 (Keck, Flat Type, Filter, Image Number)
    │       ├── KNarrowLampOnKs1.fits
    │       ├── KNarrowLampOffH1.fits  
    │       ├── KNarrowLampOffH2.fits  
    │       ├── KNarrowLampOnH1.fits 
    │       ├── KNarrowLampOnH2.fits      
    │       ├── KWideLampOnKs1.fits
    │       └── KWideLampOnKs2.fits 
    ├── Objects               
    │   └── ob151395          
    │       ├── 2018-08-05 
    │           ├── KOB151395Ks1.fits        # K ob151395 Ks 1 (Keck, Object Name, Filter, Image Number)
    │           ├── KOB151395Ks2.fits        # where Object Name is the microlensing convention:
    │           └── KOB151395Ks3.fits        # ob 15 1395 (telescope, year, event)
    │       └── 2018-08-06  
    │           ├── KOB151395NKs1.fits       # As above except the "N" denotes this is an image taken
    │           ├── KOB151395NKs2.fits       # with the NARROW camera. Images without this signifier are taken
    │           └── KOB151395NKs3.fits       # with the WIDE camera.
    │   └── mb10353          
    │       └── 2018-08-07   
    │           ├── KMB10353H1.fits          
    │           └── KMB10353H2.fits          
    ├── Skys                  
    │   ├── 2018-08-05  
    │       ├── KskyNarrowKs1.fits            # K skyNarrow Ks 1 (Keck, Sky Type, Filter, Image Number)
    │       └── KskyNarrowKs2.fits       
    ├── Source                                # Backup of Original KOA Download
    ├── calib.py
    └── tools.py

**Processed Data Folder Structure**

    ├── Darks 
    │   ├── 2018-08-05                
    │       ├── KDark30sec1.fits           
    │       ├── ...
    │       └── KDark30sec10.fits 
    │   └── KDark30sec0805.fits              # Combined Dark frame with Date.                   
    ├── Flats               
    │   ├── 2018-08-05                
    │       ├── KNarrowLampOffKs1.fits        
    │       ├── ...
    │       └── KWideLampOnKs2.fits 
    │   ├── KFlatKs0805.fits                 # Combined Flats with Date. Images are for the wide camera unless
    │   ├── KFlatH0805.fits                  # signified with an 'N'            
    │   └── KFlatNKs0805.fits                      
    ├── Objects               
    │   └── ob151395          
    │       ├── 2018-08-05 
    │           ├── Source
    │               ├── KOB151395Ks1.fits        
    │               ├── ...        
    │               ├── KOB151395Ks10.fits 
    │               ├── cor_KOB151395Ks1.fits          # STEP 1: Individual Frames Corrected for Distortion
    │               ├── ...              
    │               ├── cor_KOB151395Ks10.fits   
    │               ├── ast_cor_KOB151395Ks1.fits      # STEP 2: Astrometry Corrected Frames
    │               ├── ...              
    │               └── ast_cor_KOB151395Ks10.fits   
    │           ├── coadd_ac_KOB151395Ks.fits         # STEP 3: Coadded master frame with SWARP
    │           └── sex_cac_KOB151395Ks.tst            # STEP 4: Photometry TST File Obtained with SEXTRACTOR
    ├── Skys                  
    │   ├── 2018-08-05  
    │       ├── KskyNarrowKs1.fits 
    │       ├── ...           
    │       └── KskyNarrowKs10.fits       
    │   ├── KSkyKs0805.fits                 # Combined SkyFlats with Date. Images are for the wide camera unless
    │   ├── KSkyH0805.fits                  # signified with an 'N'
    │   └── KSkyNKs0805.fits  
    ├── Source                              
    ├── calib.py
    └── tools.py