#from pylab import *
import os
import time
import shutil
import argparse
import numpy as np
from py3shape.image import Image
from py3shape.options import Options
from py3shape.analyze import analyze, count_varied_params
from py3shape import lib
from py3shape.output import TextOutput
import sys



#I think there is no harm in doing this here, although we 
# might want to do it again later.
lib.i3_gsl_init_rng()

# Set up and parse the command line arguments using the nice Python argparse module
description = 'Im3shape measures the shapes of galaxies in astronomical survey images,\n \
taking into account that they have been distorted by a point-spread function.\n \
For more info visit https://bitbucket.org/joezuntz/im3shape\n\n \
Standard usage:\n \
parameter_filename image_file object_catalogue_file psf_file output_filename\n \
[first_image_to_process] [last_image_to_process] [additional ini file options]'

parser = argparse.ArgumentParser(description=description, add_help=True,
                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('ini_filename', type=str, help='parameter_filename')
parser.add_argument('fit_filename', type=str, help='image_file')
parser.add_argument('cat_filename', type=str, help='object_catalogue_file')
parser.add_argument('psf_filename', type=str, help='psf_file')
parser.add_argument('out_filename', type=str, help='output_filename')
parser.add_argument('first', type=int, default=0, help='[first_image_to_process]')
parser.add_argument('last', type=int, default=1e9, help='[last_image_to_process]')
#parser.add_argument('more_options', type=str, help='[additional ini file options]')



def main_203(argv):
    from py3shape.pyfits import pyfits
    args = parser.parse_args(argv)

    # Read in option file
    options = Options(args.ini_filename)
    options.validate()

    # Read in the FITS data 
    data = pyfits.getdata(args.fit_filename)

    # Read in object catalog
    gal_cat = np.atleast_2d(np.loadtxt(args.cat_filename))

    if options.psf_input == 'moffat_catalog':
        # Read in psf catalog
        psf_cat = np.loadtxt(args.psf_filename)
    elif options.psf_input == 'psf_image_single':
        # Read in psf fits image
        psf_array = pyfits.getdata(args.psf_filename)
        psf = Image(psf_array)
    elif options.psf_input == 'psf_image_cube':
        psf_arrays = pyfits.getdata(args.psf_filename)

    elif options.psf_input == "no_psf":
        psf = None
    else:
        #error!
        raise IOError('Unknown PSF input')

    ##########################################
    #   EDITED SEMENTATION MASK
    ##########################################

    # Read in segementation mask
    if options.use_segmentation_mask:
        seg_mask = pyfits.getdata(options.segmentation_mask_filename)
    else:
        seg_mask = None
    
    ##########################################

    # Create i3_image of certain stamp size
    stamp_size = options.stamp_size

    #overwrite the output filename
    #options.output_filename = args.out_filename
    options.save_output = False
    output = TextOutput(args.out_filename, args.out_filename+".epoch")
    extra_cols = []
    if options.psf_input == 'moffat_catalog':
        extra_cols.extend(('psf_beta', 'psf_e1', 'psf_e2'))
    extra_lines = [
        'driver: py3shape-main', 
        'ini: %s' % (args.ini_filename,),
        'cat: %s' % (args.cat_filename,),
        'fit: %s' % (args.fit_filename,),
        'psf: %s' % (args.psf_filename,),
        'out: %s' % (args.out_filename,),
        'first: %s' % (args.first,),
        'last: %s' % (args.last,),

        ]

    # Loop over objects in catalog
    print 'Analyzing %d/%d galaxies' % (args.last-args.first+1, len(gal_cat))
    start_time = time.clock()
    data_ymin = 0
    data_xmin = 0
    data_xmax = data.shape[1]
    data_ymax = data.shape[0]

    #main galaxy loop
    for i in range(args.first, args.last+1):
        # Read galaxy catalog entry
        identifier = gal_cat[i, 0].astype('int')
        xpos = gal_cat[i, 1]
        ypos = gal_cat[i, 2]

        # Cut out stamp from fits image
        half=stamp_size/2.

        over_edge = (ypos-half<data_ymin) | (xpos-half<data_xmin) | (ypos+half>data_ymax) | (xpos+half>data_xmax)
        if over_edge:
            print "Galaxy straddling edge of the image will not be analyzed: ", identifier
            continue
        
        stamp_array=data[int(ypos-half):int(ypos+half), int(xpos-half):int(xpos+half)]
        galaxy = Image(stamp_array)
        
        ##########################################
        #   EDITED SEMENTATION MASK
        ##########################################
        if seg_mask:
            stamp_seg_mask=seg_mask[int(ypos-half):int(ypos+half), int(xpos-half):int(xpos+half)]
            stamp_mask = np.where(stamp_seg_mask==identifier, 1, 0)
            mask = Image(stamp_mask)
        else:
            mask = None
		##########################################
		
        extra_output = {}
        if options.psf_input == 'moffat_catalog':
            # Read PSF catalog entry
            beta = psf_cat[i, 0]
            fwhm = psf_cat[i, 1]
            e1 = psf_cat[i, 2]
            e2 = psf_cat[i, 3]
            extra_output['psf_e1'] = e1
            extra_output['psf_e2'] = e2
            extra_output['psf_fwhm'] = fwhm
            extra_output['psf_beta'] = beta

            # Create PSF image from Moffat parameters
            psf = Image.make_great10_psf(stamp_size, stamp_size, beta, fwhm, e1, e2, options)
        elif options.psf_input == "psf_image_cube":
            psf = Image(psf_arrays[i])
        result, best_fit = analyze(galaxy, psf, options, mask=mask, ID=identifier, x=xpos, y=ypos)

        nparam = count_varied_params(options)        
        main, epoch = result.as_dict(0, nparam)
        main.update(extra_output)
        output.write_row(main, epoch)
        # Compute rgpp/rp
        # rgpp_rp = compute_rgpp_rp(result, psf, options)
        
        # Write out results within python rather than from the C side
        
        

    total_time = time.clock() - start_time
    ngal = args.last+1 - args.first
    print "Total time for processing in Python:" , total_time
    print "Time per galaxy:" , total_time / ngal

        

if __name__ == "__main_203__":
    main_203(sys.argv)
    
