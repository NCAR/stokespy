import numpy as np
import os
import natsort  # Sorting package
import astropy.wcs
import astropy.time
import datetime
import sunpy.map
from sunpy.net import Fido, attrs
import stokespy

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
            use_fnames.append(os.path.join(dir_path, file_))
        
    if show:
        for i, file_ in enumerate(use_fnames):
            print(i, file_)
    
    return use_fnames

def load_HMI_stokes(user_dir, user_date, user_email, max_conn=1, download=False, show_files=False, derotate=False):
    """
    Locate and fetch the HMI 720s Level 1 Stokes data closest in time to a user specified date.
    Routine can both fetch the data from JSOC or read pre-downloaded data.
    The original resolution and orientation of the data is preserved unless the derotate 
    keyword is specified which corrects for HMI's orientation relative to the pre-defined solar north direction.
    
    TODO:
    1. If download is False the code finds the nearest files in time which may not be inside the time search window.
    Create a keyword that determines the width of the search time window and then excludes any files, even locally from it.
    The downside to this is the reduced ability to load any files.
    
    Parameters
    ----------
    user_date: `astropy.time` object.
    user_email: Notification email. This must be registered with JSOC.
    user_dir: Directory where data is/will be stored.
    max_conn: The number of connections to be used when downloading from JSOC. The default setting will be the slowest but least likely to generate download errors.
    download: Flag that can be set to avoid quering JSOC for data that the user already knows is present.
    show_files: Flag to display the files where data is loaded from
    derotate: Flag if the HMI data should be derotated based on CROTA2 keyword. This ensures the PCij matrix is diagonal. 
    """
    
    # Calculate a 1s bounding time around the input date user_date
    # FIDO finds all series where at least one observation was present in the 
    # time interval.
    time0 = astropy.time.Time(user_date.gps - 1., format='gps', scale='tai')
    time1 = astropy.time.Time(user_date.gps + 1., format='gps', scale='tai')

    a_time = attrs.Time(time0, time1)
    print('Time window used for the search: ', a_time)
    
    # Set the notification email. This must be registered with JSOC. 
    a_notify = attrs.jsoc.Notify(user_email)
    
    # Set the default data directory if no user directory is specified.
    # This needs to change.
    #if user_dir is None:
    #    # Set working directory. 
    #    user_dir = os.getcwd() + '/Data/SDO/'
    #    print('User directory pointing to SDO data is not included.')
    #    print('Setting the default directory to: ' + user_dir)
    print('User data directory: ' + user_dir)
    
    # Check if the data directory exists and create one if it doesn't.
    #if not os.path.exists(user_dir):
    #    print('Data directory created: ', user_dir)
    #    os.makedirs(user_dir)
    
    ### Get the 720s HMI Stokes image series ###
    a_series = attrs.jsoc.Series('hmi.S_720s')
    
    if download:
        print('Requested data for download from JSOC.')
        results_stokes = Fido.search(a_time, a_series, a_notify)
        down_files = Fido.fetch(results_stokes, path=user_dir, max_conn=1)
        # Sort the input filenames
        all_fnames_stokes = natsort.natsorted(down_files)
    else:
        print('No download requested.')
        all_fnames_stokes = parse_folder(dir_path=user_dir, inst='hmi', series='S_720s', ext='fits', show=show_files)
        
    if len(all_fnames_stokes) > 1:
        tstamps = [i.split('.')[2] for i in all_fnames_stokes]
        tstamps = [sunpy.time.parse_time('_'.join(i.split('_')[0:2])) for i in tstamps]
        tstamps_diff = [np.abs(i.gps - user_date.gps) for i in tstamps]

    # Search for the closest timestamp
    tstamps_diff = np.asarray(tstamps_diff)
    tstamps_ix, = np.where(tstamps_diff == tstamps_diff.min())

    all_fnames_stokes = np.asarray(all_fnames_stokes)[tstamps_ix]

    # Check for individual timestamps.
    unique_tstamps = []
    for i in range(len(all_fnames_stokes)):
        if all_fnames_stokes[i].split('.')[2] not in unique_tstamps:
            unique_tstamps.append(all_fnames_stokes[i].split('.')[2])
    
    print(f'Loaded {len(all_fnames_stokes)} Stokes files with unique timestamp(s):')
    print(unique_tstamps)
    
    ## Create data array
    ## Use sunpy.map.Map to read HMI files since it provides the correct observer frame of reference.

    if derotate:
        print('OBS: Derotating each image')
    
    lvl1_data = []
    for i, fname in enumerate(all_fnames_stokes):
        tmp_map = sunpy.map.Map(fname)
        if derotate:
            tmp_map = tmp_map.rotate(order=3)
        lvl1_data.append(tmp_map.data)

    lvl1_data = np.asarray(lvl1_data)
    lvl1_data = lvl1_data.reshape(4,6,lvl1_data.shape[1], lvl1_data.shape[2])

    ## Create the WCS object
    # Expand the coordinate axis to include wavelength and stokes dimensions.

    l0 = 6173.345 * 1.e-10  # m Central wavelength for FeI line
    dl = 0.0688   * 1.e-10  # m 

    # Generate WCS for data cube using same WCS celestial information from AIA map.
    # This reads the header from the last tmp_map created (and maybe rotated) above.
    wcs_header = tmp_map.wcs.to_header()
        
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

    lvl1_wcs = astropy.wcs.WCS(wcs_header)
    
    meta = {'inst':'SDO/HMI',\
            'tobs':unique_tstamps[0],\
            'user_dir':user_dir
           }
    
    lvl1_c_HMI = stokespy.StokesCube(lvl1_data, lvl1_wcs, meta=meta)
    print(f'Created Stokes data cube with dimensions: {lvl1_c_HMI.data.shape}')
    
    return lvl1_c_HMI
    
    
def load_HMI_magvec(user_dir, user_date, user_email, max_conn=1, download=False, show_files=False, derotate=False):
    """
    Locate and fetch the HMI 720s Level 2 inversion data closest in time to a user specified date.
    Routine can both fetch the data from JSOC or read pre-downloaded data.
    The original resolution and orientation of the data is preserved unless the derotate 
    keyword is specified which corrects for HMI's orientation relative to the pre-defined solar north direction.
    
    TODO:
    1. If download == False the code finds the nearest files in time which may not be inside the time search window.
    Create a keyword that determines the width of the search time window and then excludes any files, even locally from it.
    The downside to this is the reduced ability to load any files.
    2. Clean up the outputs showing the time window better
    
    Parameters
    ----------
    user_date: `astropy.time` object.
    user_email: Notification email. This must be registered with JSOC.
    user_dir: Directory where data is/will be stored.
    max_conn: The number of connections to be used when downloading from JSOC. The default setting will be the slowest but least likely to generate download errors.
    download: Flag that can be set to avoid quering JSOC for data that the user already knows is present.
    show_files: Flag to display the files where data is loaded from
    derotate: Flag if the HMI data should be derotated based on CROTA2 keyword. This ensures the PCij matrix is diagonal. 
    """
    
    # Calculate a 1s bounding time around the input date user_date
    # FIDO finds all series where at least one observation was present in the 
    # time interval.
    time0 = astropy.time.Time(user_date.gps - 1., format='gps', scale='tai')
    time1 = astropy.time.Time(user_date.gps + 1., format='gps', scale='tai')

    a_time = attrs.Time(time0, time1)
    print('Time window used for the search: ', a_time)
    
    # Set the notification email. This must be registered with JSOC. 
    a_notify = attrs.jsoc.Notify(user_email)
    
    # Set the default data directory if no user directory is specified.
    # This needs to change.
    #if user_dir is None:
        # Set working directory. 
    #    user_dir = os.getcwd() + '/Data/SDO/'
    #    print('User directory pointing to SDO data is not included.')
    #    print('Setting the default directory to: ' + user_dir)
    print('User data directory: ' + user_dir)
    
    # Check if the data directory exists and create one if it doesn't.
    #if not os.path.exists(user_dir):
    #    print('Data directory created: ', user_dir)
    #    os.makedirs(user_dir)
    
    ### Get the HMI Milne-Eddington magentic field inversion series ###
    a_series = attrs.jsoc.Series('hmi.ME_720s_fd10')
    
    if download:
        print('Requested data for download from JSOC.')
        results_magvec = Fido.search(a_time, a_series, a_notify)
        down_files = Fido.fetch(results_magvec, path=user_dir, max_conn=1)
        # Sort the input names
        all_fnames_magvec = natsort.natsorted(down_files)
    else:
        print('No download requested.')
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

    unique_tstamps = []
    for i in range(len(all_fnames_magvec)):
        if all_fnames_magvec[i].split('.')[2] not in unique_tstamps:
            unique_tstamps.append(all_fnames_magvec[i].split('.')[2])
    
    print(f'Loaded {len(all_fnames_magvec)} inversion files with timestamps: ')
    print(unique_tstamps)
    
    # Create MagVectorCube from HMI inversions
    mag_params = ['field', 'inclination', 'azimuth']
    #mag_params = ['field']
    
    lvl2_data = []

    # Load 2D maps into lvl2_data in the order determined by entries in mag_params
    if derotate:
        print('OBS: Derotating each magnetic field image.')
    use_fnames = []
    for mag_param in mag_params:
        for i, fname in enumerate(all_fnames_magvec):
            data_id = fname.split('.')[-2]
            if data_id == mag_param:
                use_fnames.append(fname)
                tmp_map = sunpy.map.Map(fname)
                if derotate:
                    tmp_map = tmp_map.rotate(order=3)
                lvl2_data.append(tmp_map.data)
                #with astropy.io.fits.open(fname) as hdulist:
                #    lvl2_data.append(hdulist[1].data)
   
    lvl2_data = np.asarray(lvl2_data)

    #print('Filenames used: ')
    #for fname in use_fnames:
    #    print(fname)
    
    # Expand the wcs coordinates to include the magnetic field parameters.
    # Generate WCS for data cube using same WCS celestial information from the sunpy.map.
    #wcs_header = sunpy.map.Map(all_fnames_stokes[0]).wcs.to_header()
    wcs_header = tmp_map.wcs.to_header()
    
    wcs_header["WCSAXES"] = 3

    # Add Magnetic field axis.
    wcs_header["CRPIX3"] = 0
    wcs_header["CDELT3"] = 1
    wcs_header["CUNIT3"] = ''
    wcs_header["CTYPE3"] = "Parameter"
    wcs_header["CRVAL3"] = 0

    lvl2_wcs = astropy.wcs.WCS(wcs_header)
    
    meta = {'inst':'SDO/HMI',\
            'tobs':unique_tstamps[0],\
            'user_dir':user_dir}
    
    # Create the HMI Cubes
    lvl2_c_HMI = stokespy.MagVectorCube(lvl2_data, lvl2_wcs, meta=meta)
    print(f'Created magnetic field data cube with dimensions: {lvl2_c_HMI.data.shape}')
    return lvl2_c_HMI
    
def load_HinodeSP_stokes(user_dir, user_date, show_files=False):
    """
    Function that loads the lvl 1 Hinode SP observations associated with the string user_date. The data for Hinode must be downloaded independently from the Hinode website.
    
    Parameters
    ----------
    user_date: string with format "yearmmdd_hhmmss" specifying the first observation in a scan sequence. 
    user_dir: Directory where the lvl 1 data is located. 
    """
    
    ### Generate the scan file list ###
    # Assume that all the files in the indicated directory are part of the scan. 
    # TODO: Validate the scan files using the header information.
    lvl1_files = []
    user_dir = os.path.join(user_dir, 'level1', user_date[0:4], user_date[4:6], user_date[6:8], 'SP3D', user_date)
    print('Loading data from: ', user_dir)
    for file in sorted(os.listdir(user_dir)):
        if file.endswith(".fits"):
            lvl1_files.append(os.path.join(user_dir, file))
    
    # Load the first file to determine the size of the data array.
    SP_lvl1 = astropy.io.fits.open(lvl1_files[0]) 
    Nx = len(lvl1_files)
    Nstokes, Ny, Nwav = SP_lvl1[0].data.shape
    lvl1_data = np.zeros((Nx, Nstokes, Ny, Nwav))
    
    xcen_a = []
    ycen_a = []
    for ix, file in enumerate(lvl1_files):
        SP_lvl1_obj = astropy.io.fits.open(file)
        lvl1_data[ix,:,:,:] = SP_lvl1_obj[0].data
        xcen_a.append(SP_lvl1_obj['Primary'].header['XCEN'])
        ycen_a.append(SP_lvl1_obj['Primary'].header['YCEN'])
    
    # Re-arrange the data to fit the shape required by 
    # stokesCube.
    lvl1_data = lvl1_data.transpose(1, 3, 2, 0) 
    
    ### Build WCS object
    head1 = SP_lvl1['PRIMARY'].header
    
    # NOTE: the data array should be in the opposite order to the WCS, 
    # as numpy arrays are row major and wcses are Cartesian (x, y) ordered.
    # data axes order: stokes, wav, y, x
    # => wcs order   : x, y, wav, stokes
    lvl1_wcs = astropy.wcs.WCS(naxis=4)
    lvl1_wcs.wcs.ctype = ["HPLN-TAN", "HPLT-TAN", "WAVE", "STOKES"]
    lvl1_wcs.wcs.cunit = ['arcsec', 'arcsec', head1['CUNIT1'], '']
    lvl1_wcs.wcs.crpix = [(Nx+1)/2, (Ny+1)/2, head1['CRPIX1'], 0]
    lvl1_wcs.wcs.crval = [np.median(xcen_a), np.median(ycen_a), head1['CRVAL1'], 0]
    lvl1_wcs.wcs.cdelt = [head1['XSCALE'], head1['YSCALE'], head1['CDELT1'], 1]
    lvl1_wcs.wcs.set()
    
    meta = {'inst':'Hinode/SP',\
            'tobs':user_date,\
            'user_dir':user_dir}
    
    lvl1_c_SP = stokespy.StokesCube(lvl1_data, lvl1_wcs, meta=meta)
    print(f'Created Stokes data cube with dimensions: {lvl1_c_SP.data.shape}')
    
    return lvl1_c_SP
    
def load_HinodeSP_magvec(user_dir, user_date, show_files=False, magnetic_params=['Field_Strength', 'Field_Inclination', 'Field_Azimuth']):
    """
    Function that loads the Hinode SP magnetic inversion results associated with the string user_date.

    Parameters
    ----------
    user_date: string with format "yearmmdd_hhmmss" specifying the first observation in a scan sequence. 
    user_dir: Directory where the lvl1 and lvl2 data is located. We assume the data is
    separated into lvl1 and lvl2 subdirectories. 
    """
    
    ### Set the default directory substructure.
    user_dir = os.path.join(user_dir, 'level2', user_date[0:4], user_date[4:6], user_date[6:8], 'SP3D', user_date)
    print('Loading data from: ', user_dir)

    ### Read the Level 2 fit data. 
    lvl2_fname = os.path.join(user_dir, user_date + '.fits')
    lvl2_SP = astropy.io.fits.open(lvl2_fname)
    
    Ny, Nx = lvl2_SP['Field_Strength'].data.shape
    # Iterate over the list of wanted magnetic parameters.
    lvl2_data = np.zeros((len(magnetic_params), Ny, Nx))
    for i, mag_par in enumerate(magnetic_params):
        lvl2_data[i,:,:] = lvl2_SP[mag_par].data
    
    ### Build WCS objects.
    head2 = lvl2_SP['PRIMARY'].header
    
    # Calculate the center coordinates from the X,Y coordinate tables.
    x_loc = lvl2_SP[38].data
    y_loc = lvl2_SP[39].data
    
    # Computer the centers of each slit.
    Ny, Nx = y_loc.shape
    ycen_a = np.zeros(Nx)
    xcen_a = np.zeros(Nx)
    for i in range(Nx):
        ycen_a[i] = np.median(y_loc[:,i]) + head2['YSCALE']/2
        xcen_a[i] = x_loc[0,i]
    
    # NOTE: the data array should be in the opposite order to the WCS, 
    # as numpy arrays are row major and wcses are Cartesian (x, y) ordered.
    # data axes order: |B|, Binc, Bazi, y, x
    # => wcs order   : x, y, Bazi, Binc, |B|
    lvl2_wcs = astropy.wcs.WCS(naxis=3)
    lvl2_wcs.wcs.ctype = ["HPLN-TAN", "HPLT-TAN", 'Parameter']
    lvl2_wcs.wcs.cunit = ['arcsec', 'arcsec', '']
    lvl2_wcs.wcs.crpix = [(head2['NAXIS1']+1)/2, (head2['NAXIS2']+1)/2, 0]
    lvl2_wcs.wcs.crval = [np.median(xcen_a), np.median(ycen_a), 0]
    lvl2_wcs.wcs.cdelt = [head2['XSCALE'], head2['YSCALE'], 1]
    lvl2_wcs.wcs.set()
    
    meta = {'inst':'Hinode/SP',\
            'tobs':user_date,\
            'user_dir':user_dir}
    
    # Create the HMI Cubes
    lvl2_c = stokespy.MagVectorCube(lvl2_data, lvl2_wcs, meta=meta)
    print(f'Created magnetic field data cube with dimensions: {lvl2_c.data.shape}')
    return lvl2_c
