#!/home/shrey/anaconda3/bin/python
import numpy as np, pandas as pd, matplotlib.pyplot as plt, h5py, argparse, logging, wget, os, time
from tqdm import tqdm

from pycbc.types.timeseries import TimeSeries
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.filter.resample import resample_to_delta_t
from pycbc.events import ThresholdCluster
from pycbc.filter import matched_filter, resample_to_delta_t
from pycbc.vetoes import power_chisq
from pycbc.filter.matchedfilter import sigma

from snrs import SNR

####################################### Parsing #####################################################################

parser = argparse.ArgumentParser(description ='Inputs')

#Input and output files
parser.add_argument('--seg-urls', required=True, type=str, metavar ='', help ='GWOSC URLs for segments')
parser.add_argument('--temp-file', required=True, type=str, metavar ='', help ='Template hdf5 file')
parser.add_argument('--output', required=True, type=str, help ='output file')

#Controlling number of segments and templaetes for testing
parser.add_argument('--seg-nums', default = 320, type=int, metavar ='', help ='number of segments', action ='store')
parser.add_argument('--temp-select', default = True, action="store_true", help ='Remove long duration templates')
parser.add_argument('--temp-nums', default = 150, type=int, metavar ='', help ='number of templates', action ='store')

#How to calculate snr
parser.add_argument('--snr-prec', default = 2, type=int, metavar ='', help ='Precission on SNR', action ='store')
parser.add_argument('--thresh', default = 4, type=int, metavar ='', help ='Threshhold', action ='store')
parser.add_argument('--window', default = 1, type=int, metavar ='', help ='Window', action ='store')
parser.add_argument('--snr-method', default = 'new', type=str, metavar ='', help ="SNR method: 'reg' or 'new'", action ='store')

#Kind of progress update
parser.add_argument("-q", "--quite", action="store_true", help="Less Updates", default=False)
parser.add_argument("-V", "--verbose", action="store_true", help="More updates with loading screens", default=False)
parser.add_argument("--time-est", action="store_true", help="Estimate total time", default=True)

args = parser.parse_args()

############################################# Control Segments ################################################################

segments = pd.read_csv(args.seg_urls)
segments = pd.read_csv(args.seg_urls).iloc[:int(args.seg_nums)] # split number of segments

def strain_segment(urlH, urlL):
    """
    Uses wget to get strain from gwosc
    returns: H, L strains for 4096 second data
    Note: Very Slow
    """
    
    [os.remove(i) for i in os.listdir() if i[-3:] == 'tmp'] #remove any temp files that stops wget
    
    ## H
    wget.download(urlH)
    old_name = urlH.split('/')[-1]
    if os.path.isfile('Hstrain.hdf5'): # way to replace previous strain segment
        os.remove('Hstrain.hdf5')
    os.rename(old_name, 'Hstrain.hdf5')
    Hstrain = resample_to_delta_t(TimeSeries(np.array(h5py.File('Hstrain.hdf5', 'r')['strain']['Strain']), delta_t= 1/4096), 1/256)
    
    ## L
    wget.download(urlL)
    old_name = urlL.split('/')[-1]
    if os.path.isfile('Lstrain.hdf5'):
        os.remove('Lstrain.hdf5')
    os.rename(old_name, 'Lstrain.hdf5')
    Lstrain = resample_to_delta_t(TimeSeries(np.array(h5py.File('Lstrain.hdf5', 'r')['strain']['Strain']), delta_t= 1/4096), 1/256)
    
    return Hstrain, Lstrain
#################################################################################################################################

############################################# Control Templates ################################################################

hd = h5py.File(args.temp_file, 'r')
ri = list(hd.attrs['ri'])
Qi = list(hd.attrs['Qi'])
M = list(hd.attrs['M'])
Ids = list(hd.attrs['id'])

class templ:
    '''
    template
    '''
    def __init__(self, Id):
        Id = int(Id)
        self.Id = Id
        self.ri = ri[Id]
        self.Qi = Qi[Id]
        self.M = M[Id]
        self.hp = TimeSeries(np.array(hd['Hp'][str(self.Id)]), delta_t=1/256)
        self.hc = TimeSeries(np.array(hd['Hc'][str(self.Id)]), delta_t=1/256)
        
        self.hp.resize(1048576)  #resize template to 4096 seconds of data with dt = 1/256
        self.hc.resize(1048576)
               
Temps = [] # Will contain all selected templates

rem = [16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 66, 67, 68, 69, 
       70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 140, 141, 142, 
       143, 144, 145, 146, 147, 148, 149]   # Templates with almost flat strain due to very long signal duration (T >> 20 s)
for i in Ids[:args.temp_nums]:
    if i not in rem or not args.temp_select:
        Temps.append(templ(i))
#################################################################################################################################

############################################# Find Triggers ################################################################

def SNR(hp, hc, strain, PSD, lf=30, hf=256):
    '''
    Calculate newsnr for plus and cross polarisations
    '''
    
    ## Plus
    snrp = abs(matched_filter(hp, strain, psd=PSD, low_frequency_cutoff=lf, high_frequency_cutoff=hf))
    snrp = np.round(snrp, decimals=args.snr_prec)                                                                   # reduced computation
    chisqp = power_chisq(hp, strain, 26, PSD, low_frequency_cutoff=lf, high_frequency_cutoff=hf, return_bins=False)/(26*2-2)
    chisqp = np.round(chisqp, decimals=1)                                                     # reduced computation (fix precssion = 1)
    
    # chisq < 1 increased snr for injections, but also raises the noise floor.
    chisqp = np.array([i if i>1 else 1 for i in chisqp]) 
    
    snrp = snrp/chisqp
    
    ## Cross
    snrc = abs(matched_filter(hc, strain, psd=PSD, low_frequency_cutoff=lf, high_frequency_cutoff=hf))
    t = snrc.start_time
    snrc = np.round(snrc, decimals=args.snr_prec)
    chisqc = power_chisq(hc, strain, 26, PSD, low_frequency_cutoff=lf, high_frequency_cutoff=hf, return_bins=False)/(26*2-2)
    chisqc = np.round(chisqc, decimals=1)
    chisqc = np.array([i if i>1 else 1 for i in chisqc])
    snrc = snrp/chisqc
    
    ## Combine plus and cross
    snrp = TimeSeries(np.round(np.sqrt(snrp**2 + snrc**2), decimals=args.snr_prec), delta_t=1/256)
    snrp.start_time = t
    
    # return the clustered snr
    return get_prop_trig(snrp)

    
def get_prop_trig(snr, threshold=args.thresh, window=args.window):
    '''
    snr clustering
    returns: trigger times, snr
    '''
    times = snr.sample_times
    thresh = ThresholdCluster(snr)
    snrv, idx = thresh.threshold_and_cluster(window=window, threshold=threshold)
    trigs = times[idx]
    return trigs, abs(snrv)

tof = 0.0105 + 1/256*2
def find_coincs(htrigs, ltrigs):
    '''
    find coincident triggers
    faster version of pycbc_coinc_findtrigs with binary search algorithm
    returns: H, L index for coincident triggers
    '''
    
    H, L = [], []
    for i in range(len(htrigs)):
        
        ht = htrigs[i]
        low = 0
        high = len(ltrigs)-1
        keep = True

        while keep:
            mid = (low + high)//2
            if abs(ltrigs[mid]-ht) <= tof:
                keep = False
            elif high <= low:
                keep = False
                mid = -1
                if abs(ltrigs[high]-ht) <= tof:
                    mid = high       
            elif low == high - 1:
                keep = False
                if abs(ltrigs[high]-ht) <= tof:
                    mid = high
                elif abs(ltrigs[low]-ht) <= tof:
                    mid = low
                else:
                    mid = -1
            else:
                if abs(ltrigs[mid]-ht) <= tof:
                    keep = False
            
                elif ht < ltrigs[mid]:
                    high = mid
                else:
                    low = mid
        if mid != -1:
            H.append(i)
            L.append(mid)        
    return H, L
#################################################################################################################################

############################################# Run ################################################################

# output file
if os.path.isfile(args.output):
        os.remove(args.output) #replace if output file already exists
File = open(args.output, 'a')

if not args.quite:
    print("Starting Loop")
    
# option for trigger count loading screen (verbose=True)
def t_load(lst): # will return tqdm(#) or range(#) based on verbose
    return lst
if args.verbose:
    t_load = tqdm

    
lf, hf = 30, 256

for i in range(segments.shape[0]):
    
    
    t = time.time()
    seg = segments['segment'].iloc[i]                             # segment time
    text = 'Segment:' + str(seg) + '\n'                           # store result and write it into output file
    
    if not args.quite:
        print("Getting strain file for segment %s out of %d" % (i+1, segments.shape[0]))
    
    ############ strain prep ######################
    Hstrain, Lstrain = strain_segment(segments['H'].iloc[i], segments['L'].iloc[i]) # Vey Slow
    
    # PSD
    HPSD = interpolate(Hstrain.psd(4), Hstrain.delta_f)
    HPSD = inverse_spectrum_truncation(HPSD, int(4 * Hstrain.sample_rate), trunc_method='hann')
    LPSD = interpolate(Lstrain.psd(4), Lstrain.delta_f)
    LPSD = inverse_spectrum_truncation(LPSD, int(4 * Lstrain.sample_rate), trunc_method='hann')
    ##################################################
    
    if not args.quite:
        print("Starting Trigger Count")
    
    for temp, z2 in zip(Temps, t_load(range(len(Temps)))):
        
        hp, hc = temp.hp, temp.hc
        
        # get triggers
        Htrigs, Hsnr = SNR(hp, hc, Hstrain, HPSD)
        Ltrigs, Lsnr = SNR(hp, hc, Lstrain, LPSD)
        
        # get index of coincident triggers and store them in text variable
        hi, li = find_coincs(Htrigs, Ltrigs)
        Hres = [str(Hsnr[ii]) for ii in hi]
        Lres = [str(Lsnr[ii]) for ii in li]
        
        text += 'tid:' + str(temp.Id) + ', M:' + str(temp.M) + ', ri:' + str(temp.ri) + ', Qi:' + str(temp.Qi) + '\n'
        
        text += '[[' + ','.join(Hres) + ']\n' 
        text += '[' + ','.join(Lres) + ']]\n'
        
    text += '\n'
    File.write(text)
    
    if args.time_est:
        tt = time.time()-t
        print('Time Taken = %f s | Estimated remaining time = %d hrs' % (round(tt), round((segments.shape[0]-i)*tt/60/60,2)))
        
File.close()
