#!/usr/bin/env python

# Develped by Steve Caron.
# Adpated from code by Mark Vanstone:https://github.com/wmvanstone/RPiShakeCode
# Displays data for a Single RaspberryShake stations
#
# Options available:
# Display all four 4D channels
# Display phase arrival times
# Display globe with phase arrivals
#

from datetime import timedelta
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy import UTCDateTime, Stream
from obspy.geodetics import gps2dist_azimuth
from obspy.geodetics.base import locations2degrees
from matplotlib import cm
from numpy import linspace
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
import numpy as np

client = Client('http://fdsnws.raspberryshakedata.com')

# Information for the Event to Plot

LABEL = "M 7.1 - Ridgecrest, Ca "
EVT_LAT = 35.770
EVT_LON = -117.599
EVT_LATLON = (EVT_LAT, EVT_LON)
EVT_TIME = "2019-07-06 03:19:53"
EVT_Z = 8


# RaspberryShake station settings

NETWORK = 'AM'                  # AM = RaspberryShake network
STATION = "RXXXX"               # Station code of local station to plot
STA_LAT = 30.000                # + for N, - for S
STA_LON = -114.0000             # + for E, - for W
CHANNEL = 'EHZ'                 # RS Geophone Channel. 1D:SHZ - (50 sps Geophone), EHZ - (100 sps Geophone)  4D:EHZ
STA_LOCNAME = 'Anytown, Ca'



# Waveform Plot settings (Customize as needed)

DURATION = 120                 # Duration time of the waveform
TIMEOFFSET = 0                  # Time to add before start of general plot to adust placement of the waveform
ARRIVAL_TIME_ADJUST = 0     # time to add for the delay of an arriving event. Used for distant/teleseismic events
F1L = 0.1                       # Low frequency  bandpass filer
F1H = 7.0                       # High frequency bandpass filer
FIGSIZE_WAVE = 17, 8            # Figure size of the figure contaning the waveform/normal graph(s)
WAVEFORM_SIZE = 800, 400        # Individual wavefrom size
DPI = 80

# Addtional Plotting options
PLOT_4D_ACCEL = False            # If True, the the 4D acel. channels will plot with the EHZ channel
PLOT_PHASE_LINES = True         # If True, phase lines will be added to the plot
PLOT_PHASE_GLOBE = False        # If True, a globe with plotted phases will be added to the waveform. Best used with long distnce quakes.
PLOT_PS_ONLY = True             # If True, only P and S waves will be plotted. Use False for distant quakes



# Calculated constants

PHASES = ["P", "pP", "PP", "S", "Pdiff", "PKP", "PKIKP", "PcP", "ScP", "ScS",
          "PKiKP", "SKiKP", "SKP", "SKS"] # All phases. Good for distant quakes

if  PLOT_PS_ONLY:
    PHASES = ["P", "S"] # Plot only S and P good for local quakes

EVT_TIME = UTCDateTime(EVT_TIME)
STARTTIME = EVT_TIME + ARRIVAL_TIME_ADJUST
ENDTIME = STARTTIME + DURATION
STA_DIST, _, BEARING = gps2dist_azimuth(STA_LAT, STA_LON, EVT_LAT, EVT_LON)   # Station dist in m from epicenter
STA_DEGREES = locations2degrees(EVT_LAT, EVT_LON, float(STA_LAT), float(STA_LON))
CHANNEL_SET = [CHANNEL]
COLORS = [cm.plasma(x) for x in linspace(0, 0.8, len(PHASES))]                # Colors from 0.8-1.0 are not very visible
DEG_DIST_TITLE = "Km DegDist: " + str(round(STA_DEGREES, 1))

EVT_TIME = EVT_TIME.strftime('%b %d, %Y %H:%M:%S')
MODEL = 'iasp91'  # Velocity model to predict travel-times through

# Intialize Other items

FILTERLABEL = ""
PHASE_PLOT_FILENAME = ''
ACCEL_PLOT_FILENAME = ''
GLOBE_PLOT_FILENAME = ''
maxamp = []


def nospaces(string):
    out = ""
    for n in string.upper():
        if n in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789":
            out += n
        else:
            out += "_"
    return out

def pad(num):
    out = str(num)
    while len(out) < 3:
        out = "0" + out
    return out

def convert(seconds):
    seconds = seconds % (24 * 3600)
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return "%02dm:%02ds" % (minutes, seconds)

#Get the waveform data

waveform = Stream() # set up a blank stream variable

if PLOT_4D_ACCEL:
    CHANNEL_SET = ["EHZ", "ENZ", "ENN", "ENE"]

for k, channel in enumerate(CHANNEL_SET):
    try:
        st = client.get_waveforms('AM', STATION, '00', channel, starttime=STARTTIME-TIMEOFFSET, endtime=ENDTIME)
        st.merge(method=0, fill_value='latest')
        st.detrend(type='demean')
        st.filter('bandpass', freqmin=F1L, freqmax=F1H)
        FILTERLABEL = "Bandpass Filter: freqmin=" + str(F1L) + ", freqmax=" + str(F1H)
        maxvalue = (st.max())
        maxvalue = format(int(maxvalue[0]), ',d')
        maxamp.append(maxvalue)
        waveform += st
    except:
        print(STATION, channel, "off line")

# Set the waveform stats info that displayes on wavefrom plot.
# The configeration below will display as RXXXX (Channel).Location Name.Brg: XX. XXmi : xxkm  MaxAmp :xxxx

for k, wave in enumerate(waveform):
    waveform[k].stats.network = STATION + ' (' + CHANNEL_SET[k] +')'
    waveform[k].stats.station = STA_LOCNAME
    waveform[k].stats.location = " Brg: " + str(round(BEARING))+ u"\N{DEGREE SIGN} "
    waveform[k].stats.channel = str(format(int(STA_DIST/1000* 0.621371), ',d')) + "mi | " + str(format(int(STA_DIST/1000), ',d')) +"km  MaxAmp: " + maxamp[k] # Distance in km


# very big data sets may need to be decimated, otherwise the plotting routine crashes
# if DURATION * len(waveform) > 6000:
#    waveform.decimate(10, strict_length=False, no_filter=True)

# Plot the waveform

fig1 = plt.figure(figsize=(FIGSIZE_WAVE), dpi=DPI)

plt.title('Plot for '+LABEL+':  '+EVT_TIME+" UTC - lat: "+str(EVT_LAT)+" lon: "+str(EVT_LON) + " Depth: " + str(EVT_Z) + DEG_DIST_TITLE, fontsize=12, y=1.07)
plt.xticks([])
plt.yticks([])
waveform.plot(size=(WAVEFORM_SIZE), type='normal', automerge=False, equal_scale=False, fig=fig1, handle=True, starttime=STARTTIME-TIMEOFFSET, endtime=ENDTIME-TIMEOFFSET, fontsize=10, bgcolor="#F7F7F7")

if PLOT_PHASE_LINES:      # Add the arrival times and vertical lines
    plotted_arrivals = []
    for j, color in enumerate(COLORS):
        phase = PHASES[j]
        model = TauPyModel(model=MODEL)
        arrivals = model.get_travel_times(source_depth_in_km=EVT_Z, distance_in_degree=STA_DEGREES, phase_list=[phase])
        printed = False
        if arrivals:
            for i in range(len(arrivals)):
                instring = str(arrivals[i])
                phaseline = instring.split(" ")
                phasetime = float(phaseline[4])
                arrivaltime = STARTTIME - timedelta(seconds=ARRIVAL_TIME_ADJUST) + timedelta(seconds=phasetime)
                if phaseline[0] == phase and printed == False and (STARTTIME <  arrivaltime < ENDTIME):
                    axvline_arrivaltime = np.datetime64(arrivaltime)
                    legend_arrivaltime = arrivaltime.strftime("%H:%M:%S")
                    plotted_arrivals.append(tuple([round(float(phaseline[4]), 2), phaseline[0], legend_arrivaltime, color, axvline_arrivaltime]))
                    printed = True
    if plotted_arrivals:
        plotted_arrivals.sort(key=lambda tup: tup[0])   #sorts the arrivals to be plotted by arrival time
        ax = fig1.gca()
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

        if PLOT_4D_ACCEL:  #Plot the 4D accel channel wave forms
            axcheck = fig1.axes
            l = 0
            while l < (len(axcheck)-2):
                exec(   "ax" + str(l + 3) + " = plt.gcf().get_axes()[" + str(l + 1) +"]" )
                exec("trans" + str(l + 3) + " = transforms.blended_transform_factory(ax.transData, ax" + str(l + 3) + ".transAxes)")
                l += 1

        axcheck = fig1.axes
        phase_ypos = 0

        for phase_time, phase_name, arrival_legend, color_plot, arrival_vertline in plotted_arrivals:
            ax.vlines(arrival_vertline, ymin=0, ymax=1, alpha=.80, color=color_plot, ls='--', zorder=1, transform=trans)
            plt.text(arrival_vertline, .015, phase_name, c=color_plot, fontsize=11, transform=trans, horizontalalignment='center', zorder=1, bbox=dict(boxstyle="round", facecolor='#f2f2f2', ec="0.5", pad=0.24, alpha=1))
            if PLOT_4D_ACCEL:
                l = 0
                while l < (len(axcheck)-2):
                    counter = str(l + 3)
                    exec("ax" + counter + ".vlines(arrival_vertline, ymin=0, ymax=1, alpha= .80 , color=color_plot,ls='--', zorder=1, transform=trans" + counter + ")")
                   # exec("plt.text(arrival_vertline, .015, phase_name,c=color_plot, fontsize=11, transform=trans" + counter + ", horizontalalignment='center', zorder=1,bbox=dict(boxstyle='round', facecolor='#f2f2f2', ec='0.5', pad=0.24, alpha=1))")
                    l += 1
            fig_txt1 = phase_name + ' Arrival: ' + str(convert(phase_time)) +' - ' + str(arrival_legend)
            plt.figtext(.2, phase_ypos, fig_txt1, horizontalalignment='right',
                        fontsize=10, multialignment='left', c=color_plot,
                        bbox=dict(boxstyle="round", facecolor='#f2f2f2',
                                  ec="0.5", pad=0.5, alpha=1))
            phase_ypos += -.04


# Add the filter information boxes

plt.figtext(0.53, .92, FILTERLABEL, horizontalalignment='center',
            fontsize=10, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='#f2f2f2',
                      ec="0.5", pad=0.5, alpha=1))

if PLOT_PHASE_GLOBE:  # Add the inset picture of the globe at x, y, width, height, as a fraction of the parent plot
    if plotted_arrivals:
        ax1 = fig1.add_axes([0.68, 0.48, 0.4, 0.4], polar=True)
        for phase_time, phase_name, arrival_legend, color_plot, arrival_vertline in plotted_arrivals:
            arrivals = model.get_ray_paths(EVT_Z, STA_DEGREES, phase_list=[phase_name])
            if (len(arrivals)) > 1:
                m = (len(arrivals)) - 1
                arrivals = arrivals[: -m or None]  #If more than one of the same phase, use the earliers phase time
            ax1 = arrivals.plot_rays(plot_type='spherical',
                                     legend=True, label_arrivals=False,
                                     plot_all=True,
                                     show=False, ax=ax1)

        # Annotate regions of the globe
        ax1.text(0, 0, 'Solid\ninner\ncore',
                 horizontalalignment='center', verticalalignment='center',
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        ocr = (model.model.radius_of_planet -
               (model.model.s_mod.v_mod.iocb_depth +
                model.model.s_mod.v_mod.cmb_depth) / 2)
        ax1.text(np.deg2rad(180), ocr, 'Fluid outer core',
                 horizontalalignment='center',
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        mr = model.model.radius_of_planet - model.model.s_mod.v_mod.cmb_depth / 2
        ax1.text(np.deg2rad(180), mr, 'Solid mantle',
                 horizontalalignment='center',
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        rad = model.model.radius_of_planet*1.15

if PLOT_PHASE_LINES:
    PHASE_PLOT_FILENAME = '-Phases'

if PLOT_4D_ACCEL:
    ACCEL_PLOT_FILENAME = '-Accel'

if PLOT_PHASE_GLOBE:
    GLOBE_PLOT_FILENAME = '-Globe'

plt.tight_layout()

plt.savefig(nospaces(LABEL) + ACCEL_PLOT_FILENAME + PHASE_PLOT_FILENAME +
            GLOBE_PLOT_FILENAME + '.png', bbox_inches='tight')
