#!/usr/bin/env python

"""
This script is a prototype tool for looping through available instruments and instrument configurations
to calculate total system throughputs.  It uses the pandeia_spidering tools from pandeia_test to generate
a list of calculations, loops through the list using the calculation configuration to instantiate a
pandeia.engine.instrument.Instrument subclass, and then queries that instance for the total efficiency.
"""

"""
    ****  IMPORTANT ****

The files created here must match the filenames exptected in src/pandeia/ui/client/js/workbook.js

"""

import sys
import os
import numbers

import pandeia_test.pandeia_spidering as sp

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import style
style.use('ggplot')
import matplotlib.pyplot as plt
from pandeia.engine.instrument_factory import InstrumentFactory
from pandeia.engine.utils import get_key_list, get_dict_from_keys
from pandeia.engine.io_utils import read_json, write_json
from astropy.io import fits

#
#  Create all combinations of the parameters and store in calcs
#
calcs = sp.generate_calcs(telescope='jwst', allcomb=True)

# Each is a __ separated list of dictionary keys to follow
# for an entry in the calcs list. 
name_list = [
    "configuration__instrument__instrument",
    "configuration__instrument__mode",
    "configuration__instrument__aperture",
    "configuration__instrument__disperser",
    "configuration__instrument__filter"
]

#
#  Loop through the calculation combination list 
#  in order to create a figure for each.
#
for c in calcs:
    conf = c['configuration']
    iconf = conf['instrument']

    # pick out the values needed for the title and filename. skip to next if file already exists.
    # using allcomb=True covers everything and generates redundant cases that have the same throughput.
    file_list = []
    title_list = []

    # Find the instrument, mdoe, aperture, dispenser and filter based on the 
    # dictionary key hierarchy defined in "name_list"
    for name in name_list:

        # Split the name based on the __ separator
        keylist = get_key_list(name)

        # Appears to get a value following a hierarchy of 
        # dictionary keys 
        val = get_dict_from_keys(c, keylist)

        # check if it's a number and, if so, concat value with the last key for better clarity
        if isinstance(val, numbers.Number):
            val = "%s_%s" % (keylist[-1], val)

        if val is not None:
            file_list.append(val)
            title_list.append(val.upper().replace('_', '-'))
        else:
            file_list.append("none")

    # Create the plot title
    title = " ".join(title_list)

    # Create the plot output filename
    filename = './tp__' + "__".join(file_list).replace(' ', '-').replace('+', 'p').lower() + ".png"

    if os.path.isfile(filename):
        print("%s already exists. skipping..." % filename)
        continue

    # make the instrument instance, set up the wavelength vector, and get out the total efficiency.
    instrument_factory = InstrumentFactory(config=conf)

    print(iconf['instrument'], iconf['mode'], iconf['aperture'])
    try:
        wave_range = instrument_factory.get_wave_range()
    except:
        print('***  Exception with get_wave_range  ***')
        pass

    wave = np.linspace(wave_range['wmin'], wave_range['wmax'], num=500)
    eff = instrument_factory.get_total_eff(wave)

    # make the plot, save it to a file, and close things up when done.
    f = plt.figure(facecolor='white')
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize='14', colors='black')
    ax.set_ylim(0.0, 1.05)
    ax.set_axis_bgcolor('white')
    plt.ylabel("Total System Throughput", fontsize=18, color='black')
    plt.xlabel(r"$\lambda$ ($\mu$m)", fontsize=18, color='black')
    plt.plot(wave, eff, linewidth=1.5)
    plt.title(title, fontsize=20)
    plt.grid('on', which='major', color='#444444', linestyle='-')
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black') 
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')
    plt.savefig(filename)
    plt.show()
    plt.close()
