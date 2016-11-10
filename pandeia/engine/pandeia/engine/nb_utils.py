# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import numpy as np

import matplotlib
from matplotlib import style
style.use('ggplot')
matplotlib.use('nbagg')
import matplotlib.pyplot as plt

from pandeia.engine.perform_calculation import perform_calculation
from pandeia.engine.io_utils import read_json, write_json
from pandeia.engine.calc_utils import build_default_calc

from ipywidgets import widgets
from IPython.display import display, clear_output
import traitlets

def get_config(filename):
    """
    read configuration data from a JSON file and create a dict using display_strings as the keys.
    these dicts will be used to populate pull-downs with descriptive names, but provide lookup to
    the values that need to be passed back to the engine.

    Parameters
    ----------
    filename: string
        filename of JSON config file to Load

    Returns
    -------
    conf: dict
        JSON configuration data with keys swapped out for display_strings where available.
    """
    conf_data = read_json(filename)
    conf = {}
    for k in conf_data:
        if "display_string" in conf_data[k]:
            conf[conf_data[k]["display_string"]] = k
        else:
            conf[k] = k
    return conf

class WFIRST_gui(object):
    """
    create a basic GUI for WFIRST calculations using ipython widgets. to be run from within a jupyter notebook.
    """
    def __init__(self):
        self.r = {}
        self.form = widgets.VBox(width="100%", background_color="#EEE")
        self.src_form = widgets.HBox(padding='10px', width="100%", visible=False)
        self.src_select = widgets.Dropdown(description="Source type: ", options=['point', 'extended'], value='point')
        self.l1 = widgets.HTML(value="Scale length (arcsec): ", margin="5px")
        self.l2 = widgets.HTML(value="Position angle (deg): ", margin="5px")
        self.l3 = widgets.HTML(value="Ellipticity: ", margin="5px")
        self.sersic = widgets.Dropdown(description="Profile: ", options=["Gaussian", "Exponential", "de Vaucouleurs"])
        self.sersic_idx = {
            "Gaussian": 0.5,
            "Exponential": 1.0,
            "de Vaucouleurs": 4.0
        }
        self.ext_scale = widgets.BoundedFloatText(value=0.2, min=0.0, max=999999.0, width=70)
        self.ellip = widgets.BoundedFloatText(value=0.0, min=0.0, max=1.0, width=70)
        self.posang = widgets.BoundedFloatText(min=0.0, max=360.0, value=0.0, width=70)
        self.src_form.children = [self.sersic, self.l1, self.ext_scale, self.l3, self.ellip, self.l2, self.posang]
        self.src_select.on_trait_change(self.on_src_select, 'value')

        norm_form = widgets.HBox(padding='10px', width="100%")
        self.flux = widgets.BoundedFloatText(description="Source flux: ", min=0.0, max=1.0e30, value=1.0)
        self.units = widgets.Dropdown(value='ujy', options=['abmag', 'njy', 'ujy', 'mjy', 'jy'])
        atwave = widgets.HTML(value=" at ", margin='5px')
        self.wave = widgets.BoundedFloatText(min=0.1, max=99999.0, value=1.5)
        waveunits = widgets.HTML(value="microns", margin='5px')
        norm_form.children = [self.flux, self.units, atwave, self.wave, waveunits]
        self.units.on_trait_change(self.on_units_select, 'value')

        sed_form = widgets.HBox(padding='10px', width="100%")
        self.sed_select = widgets.Dropdown(
            description="SED type: ",
            options=['power-law', 'blackbody', 'star', 'extragalactic'],
            value='power-law'
        )
        self.pl_index = widgets.FloatText(description="Index: ", value=1.0, visible=True, width=50)
        self.bb_temp = widgets.BoundedFloatText(description="Temp (K): ", min=0.0, max=99999.0, value=6500.0, visible=False, width=75)

        star_config_file = os.path.join(os.environ['pandeia_refdata'], 'sed', 'phoenix', 'spectra.json')
        self.star_config = get_config(star_config_file)
        self.stars = widgets.Dropdown(options=sorted(self.star_config.keys()), visible=False)

        gal_config_file = os.path.join(os.environ['pandeia_refdata'], 'sed', 'brown', 'spectra.json')
        self.gal_config = get_config(gal_config_file)
        self.galaxies = widgets.Dropdown(options=sorted(self.gal_config.keys()), visible=False)

        self.redshift = widgets.BoundedFloatText(description="Redshift:", min=0.0, max=99999.0, value=0.0, width=70)

        self.sed_select.on_trait_change(self.on_sed_select, 'value')
        sed_form.children = [self.sed_select, self.pl_index, self.bb_temp, self.stars, self.galaxies, self.redshift]

        imager_config = read_json(os.path.join(os.environ['pandeia_refdata'], 'wfirst', 'wfirstimager', 'config.json'))
        ifu_config = read_json(os.path.join(os.environ['pandeia_refdata'], 'wfirst', 'wfirstifu', 'config.json'))

        inst_form = widgets.HBox(padding='10px', width="100%")
        self.inst_select = widgets.Dropdown(description="Instrument: ", options=['Imager'], value='Imager')
        im_filters = imager_config['filters']
        im_readmodes = imager_config['readmodes']
        im_subarrays = imager_config['subarrays']
        self.filt = widgets.Dropdown(description="Filter:", options=im_filters)
        self.readmode = widgets.Dropdown(description="Readmode:", options=im_readmodes)
        self.subarray = widgets.Dropdown(description="Sub-array:", options=im_subarrays)
        inst_form.children = [self.inst_select, self.filt, self.readmode, self.subarray]

        det_form = widgets.HBox(padding='10px', width="100%")
        self.ngroups = widgets.BoundedIntText(description="Groups: ", min=2, max=999, value=10, width=50)
        self.nints = widgets.BoundedIntText(description="Integrations: ", min=1, max=999, value=1, width=50)
        self.nexps = widgets.BoundedIntText(description="Exposures: ", min=1, max=999, value=1, width=50)
        det_form.children = [self.ngroups, self.nints, self.nexps]

        strat_form = widgets.VBox(padding='10px')
        ap_lab = widgets.HTML(value="Aperture radius (arcsec): ", margin='5px')
        self.ap_size = widgets.BoundedFloatText(min=0.0, max=999.0, value=0.1, width=60)
        self.ap_size.on_trait_change(self.check_ann, 'value')
        self.overplot = widgets.Checkbox(description="Overlay", value=True)
        self.overplot.on_trait_change(self.update_plots)
        hb1 = widgets.HBox(padding='10px', width="100%", children=[ap_lab, self.ap_size, self.overplot])
        bg_lab = widgets.HTML(value="Background annulus radii (arcsec): ", margin='5px')
        self.ann_inner = widgets.BoundedFloatText(description="inner", min=0.0, max=999.0, value=0.2, width="100%")
        self.ann_inner.on_trait_change(self.check_ann_inner, 'value')
        self.ann_outer = widgets.BoundedFloatText(description="outer", min=0.0, max=999.0, value=0.3, width="100%")
        self.ann_outer.on_trait_change(self.check_ann_outer, 'value')
        hb2 = widgets.HBox(padding='10px', width="100%", children=[bg_lab, self.ann_inner, self.ann_outer])
        strat_form.children = [hb1, hb2]

        self.calc_button = widgets.Button(description='Calculate', width="100%", background_color="#bee2c4")
        self.calc_button.on_click(self.run_calc)

        self.plot_form = widgets.HBox(padding='10px', width="100%", pack='center')
        self.oned_plots = {
            "Input Source Flux (mJy)": "target",
            "Input Background (MJy/sr)": "bg",
            "Focal Plane Rate (e-/sec/pixel)": "fp"
        }
        self.oned_units = {
            "target": "mJy",
            "bg": "MJy/sr",
            "fp": "e-/sec/pixel"
        }
        self.twod_plots = {
            "Detector (e-/sec)": "detector",
            "S/N": "snr",
            "Saturation": "saturation"
        }
        self.twod_units = {
            "detector": "e-/sec",
            "snr": "S/N",
            "saturation": ""
        }
        self.oned_pulldown = widgets.Dropdown(
            description="1D Plot",
            options=sorted(self.oned_plots.keys()),
            value="Input Source Flux (mJy)"
        )
        self.twod_pulldown = widgets.Dropdown(
            description="2D Image",
            options=sorted(self.twod_plots.keys()),
            value="S/N"
        )

        self.plot_form.children = [self.oned_pulldown, self.twod_pulldown]
        self.plot_form.visible = False
        self.oned_pulldown.on_trait_change(self.update_plots)
        self.twod_pulldown.on_trait_change(self.update_plots)

        tlab1 = widgets.HTML(value="<b>Extracted S/N: <b>", margin='5px')
        self.esn = widgets.HTML(value="0.0", margin='5px')
        tlab2 = widgets.HTML(value="       <b>Extracted Flux (e-/sec): </b>", margin='5px')
        self.eflux = widgets.HTML(value="0.0", margin='5px')
        tlab3 = widgets.HTML(value="       <b>Exposure Time (sec): <b>", margin='5px')
        self.etime = widgets.HTML(value="0.0", margin='5px')

        self.tab_form = widgets.HBox(padding='10px', width="100%", pack='center')
        self.tab_form.children = [tlab1, self.esn, tlab2, self.eflux, tlab3, self.etime]
        self.tab_form.visible = False
        self.form.children = [
            self.src_select,
            self.src_form,
            norm_form,
            sed_form,
            inst_form,
            det_form,
            strat_form,
            self.calc_button,
            self.tab_form,
            self.plot_form
        ]

    def update_plots(self):
        """
        update the 1D and 2D plots.  they're part of the same figure so have to be drawn together.  hard to do
        two independent plots in one cell in a notebook.
        """
        oned_key = self.oned_plots[self.oned_pulldown.value]
        twod_key = self.twod_plots[self.twod_pulldown.value]
        oned_curve = self.r['1d'][oned_key]
        twod_im = self.r['2d'][twod_key]
        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_subplot(121)
        plot = ax1.plot(oned_curve[0], oned_curve[1])
        ax1.set_xlabel(r'$\mu m$')
        ax1.set_ylabel(self.oned_units[oned_key])
        t = self.r['transform']
        xmin = t['x_min']
        xmax = t['x_max']
        ymin = t['y_min']
        ymax = t['y_max']
        extent = [xmin, xmax, ymin, ymax]
        ax2 = fig.add_subplot(122)
        ax2.set_xlabel("arcsec")
        ax2.set_ylabel("arcsec", rotation=270)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.set_ticks_position("both")
        im = plt.imshow(twod_im, interpolation='nearest', extent=extent)
        if self.overplot.value is True:
            circles = []
            circles.append(plt.Circle((0, 0), radius=self.ap_size.value, edgecolor='white', facecolor='none'))
            circles.append(plt.Circle((0, 0), radius=self.ann_inner.value, edgecolor='red', facecolor='none'))
            circles.append(plt.Circle((0, 0), radius=self.ann_outer.value, edgecolor='red', facecolor='none'))
            for c in circles:
                im.axes.add_artist(c)
        if twod_key == "saturation":
            norm = matplotlib.colors.Normalize(vmin=0, vmax=2)
            im.set_norm(norm)
            c = plt.colorbar(ax=ax2, orientation="horizontal", label=self.twod_units[twod_key], ticks=[0, 1, 2])
            c.ax.set_xticklabels(["None", "Soft", "Hard"])
        else:
            c = plt.colorbar(ax=ax2, orientation="horizontal", label=self.twod_units[twod_key])
        plt.tight_layout()
        plt.show()
        clr = clear_output(wait=True)

    @property
    def display(self):
        """
        display the GUI
        """
        display(self.form)

    @property
    def calc_results(self):
        """
        return the calculation results
        """
        return self.r

    @calc_results.setter
    def calc_results(self, r):
        self.r = r
        self.update_plots()

    def check_ann(self, name, value):
        """
        check the background estimation annulus to make sure it's valid

        Parameters
        ----------
        name: string
            not used, but expected for on_trait_change callbacks
        """
        self.check_ann_inner(name=name, value=value)
        self.check_ann_outer(name=name, value=value)

    def check_ann_inner(self, name, value):
        if self.ann_inner.value <= self.ap_size.value:
            self.ann_inner.value = round(self.ap_size.value + 0.1, 3)
        self.check_ann_outer(name=name, value=value)

    def check_ann_outer(self, name, value):
        if self.ann_outer.value - self.ann_inner.value <= 0.1:
            self.ann_outer.value = round(self.ann_inner.value + 0.1, 3)

    def on_src_select(self, name, value):
        if value == 'point':
            self.src_form.visible = False
        else:
            self.src_form.visible = True

    def on_units_select(self, name, value):
        if value == 'abmag':
            self.flux.value = 25.0
        else:
            self.flux.value = 1.0

    def on_sed_select(self, name, value):
        if value == 'power-law':
            self.pl_index.visible = True
            self.bb_temp.visible = False
            self.stars.visible = False
            self.galaxies.visible = False
        elif value == 'blackbody':
            self.pl_index.visible = False
            self.bb_temp.visible = True
            self.stars.visible = False
            self.galaxies.visible = False
        elif value == 'star':
            self.pl_index.visible = False
            self.bb_temp.visible = False
            self.stars.visible = True
            self.galaxies.visible = False
        elif value == 'extragalactic':
            self.pl_index.visible = False
            self.bb_temp.visible = False
            self.stars.visible = False
            self.galaxies.visible = True
        else:
            self.pl_index.visible = False
            self.bb_temp.visible = False
            self.stars.visible = False
            self.galaxies.visible = False

    def run_calc(self, b):
        c = build_default_calc("wfirst", "wfirstimager", "imaging")
        c['configuration']['detector']['nexp'] = self.nexps.value
        c['configuration']['detector']['ngroup'] = self.ngroups.value
        c['configuration']['detector']['nint'] = self.nints.value
        c['configuration']['detector']['readmode'] = self.readmode.value
        c['configuration']['detector']['subarray'] = self.subarray.value
        c['configuration']['instrument']['filter'] = self.filt.value

        src = c['scene'][0]
        if self.src_select.value == "extended":
            src['shape']['geometry'] = 'sersic'
            a = self.ext_scale.value
            e = self.ellip.value
            b = (1.0 - e) * a
            s_idx = self.sersic_idx[self.sersic.value]
            # if gaussian, convert a/b to sigma
            if s_idx == 0.5:
                a *= np.sqrt(2.0)
                b *= np.sqrt(2.0)
            src['shape']['major'] = a
            src['shape']['minor'] = b
            src['shape']['sersic_index'] = s_idx
            src['position']['orientation'] = self.posang.value

        src['spectrum']['redshift'] = self.redshift.value
        src['spectrum']['normalization']['norm_flux'] = self.flux.value
        src['spectrum']['normalization']['norm_fluxunit'] = self.units.value
        src['spectrum']['normalization']['norm_wave'] = self.wave.value

        sed = self.sed_select.value
        if sed == "power-law":
            src['spectrum']['sed']['sed_type'] = "powerlaw"
            src['spectrum']['sed']['index'] = self.pl_index.value
        if sed == "blackbody":
            src['spectrum']['sed']['sed_type'] = "blackbody"
            src['spectrum']['sed']['temp'] = self.bb_temp.value
        if sed == "star":
            src['spectrum']['sed']['sed_type'] = "phoenix"
            src['spectrum']['sed']['key'] = self.star_config[self.stars.value]
        if sed == "extragalactic":
            src['spectrum']['sed']['sed_type'] = "brown"
            src['spectrum']['sed']['key'] = self.gal_config[self.galaxies.value]

        c['strategy']['aperture_size'] = self.ap_size.value
        ann = [self.ann_inner.value, self.ann_outer.value]
        c['strategy']['sky_annulus'] = ann

        self.r = perform_calculation(c, dict_report=True)
        self.calc_input = c
        self.plot_form.visible = True
        self.esn.value = "%.2f" % self.r['scalar']['sn']
        self.eflux.value = "%.2f" % self.r['scalar']['flux']
        self.etime.value = "%.2f" % self.r['scalar']['on_source_time']
        self.tab_form.visible = True

        self.update_plots()
