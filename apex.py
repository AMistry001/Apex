"""
Automatic Peak Extractor - APEx

Coded by: Adam Greenberg - adamhgreenberg@astro.ucla.edu
"""
import numpy as np
import argparse as ap
import sys

from gooey import Gooey

#REMEMBER TO REMOVE +1 FROM trim_right AND TO UNDO CHANGES MADE TO READTIFF IN PUTIL

# note that a profile function f must always be of the form:
# 	f(xs, height, loc, width)
profiles = {
		"none"	: None,
		"gauss" : lambda xs,height,loc,width: height*np.exp( -0.5*( (xs-loc)/width )**2),
}
@Gooey()
def main():

    parser = ap.ArgumentParser(
    		description="extracts peak locations and ancillary information from an unrolled diffraction file.",
    		formatter_class=ap.ArgumentDefaultsHelpFormatter
    		)
    
    parser.add_argument("-s","--std_threshold",type = float, default = 2.0, help = "minimum threshold for peak detection, in units of standard deviation of the noise")
    parser.add_argument("-t","--trim_left",type = int, default = 960, help = "number of columns to discard from left of image")
    parser.add_argument("-T","--trim_right",type = int, default = 2048 - 128, help = "number of columns to discard from right of image")
    parser.add_argument("-w","--window_size",type = int, default = 32, help = "window size (in pixels) for background detection/interpolation")
    parser.add_argument("-n","--neighbors",type = int, default = 4, help = "minimum number of smaller neighbors needed for peak detection")
    parser.add_argument("-a","--athresh",type = int, default = 4.0, help = "scale threshold for wavelet smoothing")
    parser.add_argument("-c","--chunk_size",type = int, default = 4, help = "number of (summed) rows per chunk")
    parser.add_argument("-f","--peak_func",type = lambda x: x.lower(), default = "none", help = "function for profile of peak. current options are: %s."%",".join(profiles.keys()))
    parser.add_argument("-N","--num_waves",type = int, default = 32, help = "number of wavelet scales to use")
    
    
    parser.add_argument("-p","--prefix",type = str, default = "", help = "prefix for output files. default uses input filename root")
    
    parser.add_argument("--verbose",default=False,action="store_true", help="ouput progress")
    parser.add_argument("--show",default=False,action="store_true", help="show plots interactively with pyplot")
    parser.add_argument("--blink",default=False,action="store_true", help="make an overlay gif")
    
    
    parser.add_argument("input_file",type = str, help = "unrolled  file")
    
    args = parser.parse_args()
    
    import matplotlib
    if not args.show: matplotlib.use("agg")
    import matplotlib.pyplot as plt
    import putil as pu
    
    if not args.prefix:
    	args.prefix = args.input_file[ : args.input_file.rfind(".")]
    
    image,chunks,peaklists,widthlists,heightlists = pu.processImage(	filename = args.input_file,
    									peakfunc = profiles[args.peak_func],
    									ltrim = args.trim_left,
    									rtrim = args.trim_right,
    									corrlen = args.chunk_size,
    									stdthresh = args.std_threshold,
    									bgwindowsize = args.window_size,
    									athresh = args.athresh,
    									numneighbors = args.neighbors,
    									nwaves = args.num_waves,
    									verbose = args.verbose, )
    
    
    eps = 0.02
    imextent = [args.trim_left,args.trim_left + image.shape[1],image.shape[0],0]
    
    figh = plt.figure(figsize = (10,10))
    figw = plt.figure(figsize = (10,10))
    axh = figh.add_subplot(111)
    axw = figw.add_subplot(111)
    if args.verbose: print("writing output data file")
    outfile = open("%s_list.txt"%args.prefix,"w")
    outfile.write("#row_start,row_end,height,peakloc,width\n")
    for chunk,heightlist,peaklist,widthlist in zip(chunks,heightlists,peaklists,widthlists):
    	axh.scatter([peak+args.trim_left for peak in peaklist],heightlist)
    	axw.scatter([peak+args.trim_left for peak in peaklist],widthlist)
    	for height,peak,width in zip(heightlist,peaklist,widthlist):
    		y0,y1=chunk
    		peak+=args.trim_left
    		outfile.write("%5i,%5i,% 10.2f,% 7.2f,% 7.2f\n"%(y0,y1,height,peak,width))
    outfile.close()
    
    if args.verbose: print("plotting heights")
    axh.set_xlim(imextent[0],imextent[1])
    axh.set_ylim(ymin=1)
    axh.set_yscale('log')
    axh.set_ylabel("peak height")
    axh.set_xlabel("column (pixels)")
    figh.savefig("%s_heights.png"%args.prefix)
    
    if args.verbose: print("plotting widths")
    axw.set_xlim(imextent[0],imextent[1])
    axw.set_ylim(ymin=0)
    axw.set_ylabel("peak width")
    axw.set_xlabel("column (pixels)")
    figw.savefig("%s_widths.png"%args.prefix)
    if args.show: plt.show()
    
    if args.verbose: print("plotting raw image")
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.imshow(np.log(image), interpolation="none", aspect="auto", cmap=plt.get_cmap("Spectral"), extent=imextent)
    plt.tight_layout(rect=[eps,eps,1-eps,1-eps])
    plt.xlim(imextent[:2])
    plt.ylim(imextent[2:][::1]) #change 1 to -1 to flip overlay image
    plt.savefig("%s_raw.png"%args.prefix)
    if args.verbose: print("plotting overlay")
    for chunk,peaklist,widthlist in zip(chunks,peaklists,widthlists):
    	for peak,width in zip(peaklist,widthlist):
    		y0,y1 = chunk
    		peak+=args.trim_left
    		ax.fill_between((peak-width/2.,peak+width/2.),(y0,y0),(y1,y1), color="b",alpha=0.3)
    plt.savefig("%s_overlay.png"%args.prefix)
    if args.show: plt.show()
    plt.clf()
    
    if args.verbose: print("plotting collapse")
    M = pu.getImageMask(image.shape, chunks,peaklists,widthlists)
    collapse = np.sum(M,axis=0)
    plt.plot(range(*imextent[:2]),collapse)
    plt.xlim(imextent[0],imextent[1])
    plt.tight_layout(rect=[eps,eps,1-eps,1-eps])
    plt.xlabel("column (pixels)")
    plt.ylabel("Number of peak detections")
    plt.savefig("%s_collapse.png"%args.prefix)
    if args.show: plt.show()
    plt.clf()
    
    if args.blink:
    	import subprocess as sp
    	cmd = "convert -delay 100 %s_overlay.png %s_raw.png %s.gif"%tuple(3*[args.prefix])
    	sp.Popen(cmd,shell=True)
        
   
if __name__ == "__main__":
    main()