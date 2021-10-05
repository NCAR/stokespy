import numpy as np
import os
import natsort  # Sorting package
import astropy.wcs
import astropy.time
import datetime
import sunpy.map
from sunpy.net import Fido, attrs

def parse_folder(dir_path=None, inst=None, wave=None, ext=None, 
                 series=None, repo='JSOC', show=True):
    """
    Search all the filenames in a folder containing SDO data and use the keywords
    to select a desired subset. For example setting inst="hmi" will obtain all the hmi filenames.
    dir_path: Path to the directory containing the data.
    inst: SDO instrument
    wave: wavelength (primarily for AIA data)
    ext: Select the file extension
    series: String characterizing the data series (e.g. hmi.S_720s, aia.lev1_euv_12s)
    repo: Choose the data repository. Each repository stores filenames with different syntaxes.
    show: Flag that enumerates the files found.
    """
    
    if dir_path is None:
        dir_path = os.getcwd()
    
    # Read and sort all the filenames in the folder.
    all_fnames = natsort.natsorted(os.listdir(dir_path)) 
    
    # Select a subset of files
    use_fnames = []
    
    for i, file_ in enumerate(all_fnames):
        if repo == 'JSOC':
            sfile_ = file_.lower().split(".")
        elif repo == 'VSO':
            sfile_ = file_.lower().split("_")
        
        if sfile_[0] == inst.lower() and sfile_[1] == series.lower(): 
            use_fnames.append(dir_path + file_)
        
        '''
        if ext is None and inst is None:
            use_fnames.append(dir_path + file_)
        elif ext is None and inst is not None and sfile_[0] == inst.lower() and wave is None:
            use_fnames.append(dir_path + file_)
        elif ext is None and inst is not None and sfile_[0] == inst.lower() and wave is not None and sfile_[2] == wave + 'a':
            use_fnames.append(dir_path + file_)
        elif ext is not None and file_.endswith(ext) and inst is None:
            use_fnames.append(dir_path + file_)
        elif ext is not None and file_.endswith(ext) and inst is not None and sfile_[0] == inst.lower() and wave is None:
            use_fnames.append(dir_path + file_)
        elif ext is not None and file_.endswith(ext) and inst is not None and sfile_[0] == inst.lower() and wave is not None and sfile_[2] == wave + 'a':
            use_fnames.append(dir_path + file_)
        '''
    '''
    for i, file_ in enumerate(all_fnames):
        # Select by extension.
        if ext is not None and file_.endswith(ext):        
            sfile_ = file_.split("_")
            # Select by instrument.
            if inst is not None and sfile_[0] == inst.lower(): 
                if wave is not None and inst.lower() == 'aia' and sfile_[2] == wave + 'a':   
                    use_fnames.append(dir_path + file_)
        else:
            sfile_ = file_.split("_")
            # Select by instrument.
            if inst is not None and sfile_[0] == inst.lower(): 
                if wave is not None and inst.lower() == 'aia' and sfile_[2] == wave + 'a':   
                    use_fnames.append(dir_path + file_)
    '''
    
    if show:
        for i, file_ in enumerate(use_fnames):
            print(i, file_)
    
    return use_fnames

def get_HMI_data(user_date, user_notify='gdima@hawaii.edu', user_dir=None, max_conn=1, download=False, show_files=False):
    """
    Locate the nearest (relative to the user input date) HMI 720s Level 1 Stokes 
    data series and Level 2 inversion results. 
    If the data is already present in the user specified data directory is loaded from there.
    
    TODO:
    1. If download=False the code finds the nearest files in time which may not be inside the time search window.
    Create a keyword that determines the width of the search time window and then excludes any files, even locally from it.
    The downside to this is the reduced ability to load any files.
    
    Parameters
    ----------
    user_date: `astropy.time` object.
    user_notify: Notification email. This must be registered with JSOC.
    
    """
    
    # Calculate a 1s bounding time around the input date user_date
    # FIDO finds all series where at least one observation was present in the 
    # time interval.
    time0 = astropy.time.Time(user_date.gps - 1., format='gps', scale='tai')
    time1 = astropy.time.Time(user_date.gps + 1., format='gps', scale='tai')

    a_time = attrs.Time(time0, time1)
    print('Time window used for the search: ', a_time)
    
    # Set the notification email. This must be registered with JSOC. 
    a_notify = attrs.jsoc.Notify(user_notify)
    
    # Set the default data directory if no user directory is specified.
    if user_dir is None:
        # Set working directory. 
        user_dir = os.getcwd() + '/Data/SDO/'
        print('User directory pointing to SDO data is not included.')
        print('Setting the default directory to: ' + user_dir)
    
    # Check if the data directory exists and create one if it doesn't.
    if not os.path.exists(user_dir):
        print('Data directory created: ', user_dir)
        os.makedirs(user_dir)
    
    ### Get the 720s HMI Stokes image series ###
    a_series = attrs.jsoc.Series('hmi.S_720s')
    
    if download:
        results_stokes = Fido.search(a_time, a_series, a_notify)
        down_files = Fido.fetch(results_stokes, path=user_dir, max_conn=1)
        # Sort the input filenames
        all_fnames_stokes = natsort.natsorted(down_files)
    else:
        all_fnames_stokes = parse_folder(dir_path=user_dir, inst='hmi', series='S_720s', ext='fits', show=show_files)
    
        if len(all_fnames_stokes) > 1:
            tstamps = [i.split('.')[2] for i in all_fnames_stokes]
            tstamps = [sunpy.time.parse_time('_'.join(i.split('_')[0:2])) for i in tstamps]
            tstamps_diff = [np.abs(i.gps - user_date.gps) for i in tstamps]
       
        # Search for the closest timestamp
        tstamps_diff = np.asarray(tstamps_diff)
        tstamps_ix, = np.where(tstamps_diff == tstamps_diff.min())

        all_fnames_stokes = np.asarray(all_fnames_stokes)[tstamps_ix]
        
        print(f'No download requested. Nearest {len(all_fnames_stokes)} Stokes files found: ')
        print(all_fnames_stokes)
    
    ### Get the HMI Milne-Eddington magentic field inversion series ###
    a_series = attrs.jsoc.Series('hmi.ME_720s_fd10')
    
    if download:
        results_magvec = Fido.search(a_time, a_series, a_notify)
        down_files = Fido.fetch(results_magvec, path=user_dir, max_conn=1)
        # Sort the input names
        all_fnames_magvec = natsort.natsorted(down_files)
    else:
        all_fnames_magvec = parse_folder(dir_path=user_dir, inst='hmi', series='ME_720s_fd10', ext='fits', show=show_files)
        
        if len(all_fnames_magvec) > 1:
            tstamps = [i.split('.')[2] for i in all_fnames_magvec]
            tstamps = [sunpy.time.parse_time('_'.join(i.split('_')[0:2])) for i in tstamps]
            tstamps_diff = [np.abs(i.gps - user_date.gps) for i in tstamps]
        else:
            print('No files found close to the date requested')
            return 
            
        # Search for the closest timestamp
        tstamps_diff = np.asarray(tstamps_diff)
        tstamps_ix, = np.where(tstamps_diff == tstamps_diff.min())

        all_fnames_magvec = np.asarray(all_fnames_magvec)[tstamps_ix]
        
        print(f'No download requested. Nearest inversion {len(all_fnames_magvec)} files found: ')
        print(all_fnames_magvec)

    ## Create data array
    ## Use sunpy.map.Map to read HMI files since it provides the correct observer frame of reference.

    level1_data = []
    for i, fname in enumerate(all_fnames_stokes):
        level1_data.append(sunpy.map.Map(fname).data)

    level1_data = np.asarray(level1_data)
    level1_data = level1_data.reshape(4,6,level1_data.shape[1], level1_data.shape[2])

    print(f'Created data cube with dimensions: {level1_data.shape}')    
    
    ## Create the WCS object
    # Expand the coordinate axis to include wavelength and stokes dimensions.

    l0 = 6173.345 * 1.e-10  # m Central wavelength for FeI line
    dl = 0.0688   * 1.e-10  # m 

    # Generate WCS for data cube using same WCS celestial information from AIA map.
    wcs_header = sunpy.map.Map(all_fnames_stokes[0]).wcs.to_header()

    wcs_header["WCSAXES"] = 4

    # Add wavelength axis.
    wcs_header["CRPIX3"] = 3.5
    wcs_header["CDELT3"] = dl
    wcs_header["CUNIT3"] = 'm'
    wcs_header["CTYPE3"] = "WAVE"
    wcs_header["CRVAL3"] = l0

    # Add Stokes axis.
    wcs_header["CRPIX4"] = 0
    wcs_header["CDELT4"] = 1
    wcs_header["CUNIT4"] = ''
    wcs_header["CTYPE4"] = "STOKES"
    wcs_header["CRVAL4"] = 0

    level1_wcs = astropy.wcs.WCS(wcs_header)
    
    # Create MagVectorCube from HMI inversions
    mag_params = ['field', 'inclination', 'azimuth']

    level2_data = []

    # Load 2D maps into level2_data in the order determined by entries in mag_params
    use_fnames = []
    for mag_param in mag_params:
        for i, fname in enumerate(all_fnames_magvec):
            data_id = fname.split('.')[-2]
            if data_id == mag_param:
                use_fnames.append(fname)
                with astropy.io.fits.open(fname) as hdulist:
                    level2_data.append(hdulist[1].data)

    level2_data = np.asarray(level2_data)

    print(f'Created data cube with dimensions: {level2_data.shape}')
    print('Filenames used: ')
    for fname in use_fnames:
        print(fname)
    
    # Expand the wcs coordinates to include the magnetic field parameters.

    # Generate WCS for data cube using same WCS celestial information from AIA map.
    wcs_header = sunpy.map.Map(all_fnames_stokes[0]).wcs.to_header()

    wcs_header["WCSAXES"] = 3

    # Add Magnetic field axis.
    wcs_header["CRPIX3"] = 0
    wcs_header["CDELT3"] = 1
    wcs_header["CUNIT3"] = ''
    wcs_header["CTYPE3"] = "Parameter"
    wcs_header["CRVAL3"] = 0

    level2_wcs = astropy.wcs.WCS(wcs_header)
    
    return level1_data, level1_wcs, level2_data, level2_wcs

def get_SP_data(user_date, user_dir=None, show_files=False):
    """
    Function that loads the Hinode SP observations associated with the the string user_date.
    ----------
    user_date: string with format "yearmmdd_hhmmss" specifying the first observation in a scan sequence. 
    user_dir: Directory where the Level1 and Level2 data is located. We assume the data is
    separated into Level1 and Level2 subdirectories. 
    """
    
    ### Set the default directory substructure.
    if user_dir is None:
        user_dir = os.getcwd() + '/Data/Hinode/'
   
    ### Generate the scan file list.
    level1_dir = user_dir + 'Level1/' + user_date
    level1_files = []
    for file in sorted(os.listdir(level1_dir)):
        if not file.endswith(".fits"):
            continue
        level1_files.append(os.path.join(level1_dir, file))
    
    SP_level1 = astropy.io.fits.open(level1_files[0])
    #SP_level1.info()

    ### Read all the Level 1 scan data 
    Nx = len(level1_files)
    Nstokes, Ny, Nwav = SP_level1[0].data.shape
    level1_data = np.zeros((Nx, Nstokes, Ny, Nwav))

    for ix, file in enumerate(level1_files):
        SP_lvl1_obj = astropy.io.fits.open(file)
        level1_data[ix] = SP_lvl1_obj[0].data

    level1_data = level1_data.transpose(1, 3, 2, 0) # data axes order: stokes, wav, y, x
    
    ### Read the Level 2 fit data. 
    level2_fname = user_dir + '/Level2/' + user_date + '.fits'
    SP_level2 = astropy.io.fits.open(level2_fname)
    
    Ny, Nx = SP_level2['Field_Strength'].data.shape
    level2_data = np.zeros((3, Ny, Nx))
    level2_data[0] = SP_level2['Field_Strength'].data
    level2_data[1] = SP_level2['Field_Inclination'].data
    level2_data[2] = SP_level2['Field_Azimuth'].data
    level2_data.shape
    
    ### Build WCS objects.
    head1 = SP_level1['PRIMARY'].header
    head2 = SP_level2['PRIMARY'].header
    
    # NOTE: the data array should be in the opposite order to the WCS, 
    # as numpy arrays are row major and wcses are Cartesian (x, y) ordered.
    # data axes order: stokes, wav, y, x
    # => wcs order   : x, y, wav, stokes
    level1_wcs = astropy.wcs.WCS(naxis=4)
    level1_wcs.wcs.ctype = ["HPLN-TAN", "HPLT-TAN", "WAVE", "STOKES"]
    level1_wcs.wcs.cunit = ['arcsec', 'arcsec', head1['CUNIT1'], '']
    level1_wcs.wcs.crpix = [(head2['NAXIS1']+1)/2, (head2['NAXIS2']+1)/2, head1['CRPIX1'], 0]
    level1_wcs.wcs.crval = [head2['XCEN'], head2['YCEN'], head1['CRVAL1'], 0]
    level1_wcs.wcs.cdelt = [head2['XSCALE'], head2['YSCALE'], head1['CDELT1'], 1]
    level1_wcs.wcs.set()
    
    # NOTE: the data array should be in the opposite order to the WCS, 
    # as numpy arrays are row major and wcses are Cartesian (x, y) ordered.
    # data axes order: |B|, Binc, Bazi, y, x
    # => wcs order   : x, y, Bazi, Binc, |B|
    level2_wcs = astropy.wcs.WCS(naxis=3)
    level2_wcs.wcs.ctype = ["HPLN-TAN", "HPLT-TAN", 'Parameter']
    level2_wcs.wcs.cunit = ['arcsec', 'arcsec', '']
    level2_wcs.wcs.crpix = [(head2['NAXIS1']+1)/2, (head2['NAXIS2']+1)/2, 0]
    level2_wcs.wcs.crval = [head2['XCEN'], head2['YCEN'], 0]
    level2_wcs.wcs.cdelt = [head2['XSCALE'], head2['YSCALE'], 1]
    level2_wcs.wcs.set()
    
    return level1_data, level1_wcs, level2_data, level2_wcs