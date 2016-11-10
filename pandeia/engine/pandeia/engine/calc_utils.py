# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division, absolute_import

import os
import pkg_resources

from . import config
from .io_utils import read_json
from pandeia.engine.source import Source, Line
from pandeia.engine.custom_exceptions import EngineInputError

default_refdata_directory = config.default_refdata_directory


def get_telescope_config(telescope):
    """
    Get the telescope configuration to tell us instruments and modes.
    Must have the pandeia_refdata env variable set and pointing to the reference data

    Parameters
    ----------
    telescope: string
        Telescope identifier string. Currently supports jwst, wfirst, and hst.

    Returns
    -------
    config: dict
        Telescope configuration data
    """
    config_file = os.path.join(default_refdata_directory, telescope, 'telescope', 'config.json')
    config = read_json(config_file, raise_except=True)
    return config


def get_instrument_config(telescope, instrument):
    """
    Get the instrument configuration.
    Must have the pandeia_refdata env variable set and pointing to the reference data

    Parameters
    ----------
    telescope: string
        Telescope identifier string. Currently supports jwst, wfirst, and hst.
    instrument: string
        Instrument identifier string.

    Returns
    -------
    config: dict
        Instrument configuration data
    """
    config_file = os.path.join(default_refdata_directory, telescope, instrument, 'config.json')
    config = read_json(config_file, raise_except=True)
    return config


def build_default_source(geometry="point"):
    """
    Build and return a default Source in engine API dict format

    Returns
    -------
    src: dict
        Source configuration data in engine API format
    """
    src = Source(shape={"geometry": geometry}).as_dict()
    return src


def build_default_line():
    """
    Build and return a default Line in engine API dict format

    Returns
    -------
    line: dict
        Line configuration data in engine API format
    """
    line = Line().as_dict()
    return line


def build_default_scene():
    """
    Build and return a default scene which consists of a single default Source

    Returns
    -------
    scene: list
        Single element list containing a source as built by build_default_source()
    """
    scene = [build_default_source()]
    return scene


def build_empty_scene():
    """
    Build and return an empty scene. Because of the way ConvolvedSceneCube's and AstroSpectrum's are created,
    there must be some concept of a spectrum and thus a source contained within a scene. The (admittedly hacky)
    workaround is to build a default scene with a flat-spectrum point source and set its flux to 0. This will provide
    the framework necessary to build things from this scene using the existing code, but not add any actual flux. The
    primary motivation for this is the BNC which needs to build a ConvolvedSceneCube that only contains background signal.

    Returns
    -------
    scene: list
        Single element list containing a zero-flux source
    """
    scene = build_default_scene()
    scene[0]['spectrum']['normalization']['norm_flux'] = 0.0
    return scene


def build_default_calc(telescope, instrument, mode, **kwargs):
    """
    Build a default calculation given an telescope, instrument, and mode

    Parameters
    ----------
    telescope: string
        Telescope identifier string. Currently supports jwst, wfirst, and hst.
    instrument: string
        Instrument identifier string.
    mode: string
        Mode identifier string.

    Returns
    -------
    calc: dict
        Engine output API compliant dict defining a calculation
    """
    calc = dict()
    ex_args = dict(**kwargs)

    ins_config = get_instrument_config(telescope=telescope, instrument=instrument)
    if mode is None or mode not in ins_config['modes']:
        mode = ins_config['defaults']['mode']

    calc['configuration'] = ins_config['defaults'][mode]

    calc['scene'] = build_default_scene()

    calc['background'] = "medium"

    st_defaults = ins_config['strategy_config'][mode]
    calc['strategy'] = dict()

    if 'method' in ex_args:
        calc['strategy']['method'] = ex_args['method']
    else:
        calc['strategy']['method'] = st_defaults['method']

    method_defaults_file = "%s.json" % calc['strategy']['method']

    method_defaults = read_json(os.path.join(default_refdata_directory, "strategy", method_defaults_file))
    for k in method_defaults:
        if 'permitted' not in k:
            calc['strategy'][k] = method_defaults[k]

    calc_defaults_file = "defaults/calculationconfig_defaults.json"
    calc['calculation'] = read_json(pkg_resources.resource_filename("pandeia.engine", calc_defaults_file))

    # sometimes aperture_key is a dict because some stratgies will use it, but others won't.
    # if it's not used, it will be set to None
    if 'aperture_key' in st_defaults:
        if isinstance(st_defaults['aperture_key'], dict):
            ap_key = st_defaults['aperture_key'][calc['strategy']['method']]
        else:
            ap_key = st_defaults['aperture_key']
        if ap_key is not None:
            key = calc['configuration']['instrument'][ap_key]
            calc['strategy']['aperture_size'] = st_defaults['aperture_sizes'][key]
    else:
        ap_key = None
        if 'aperture_size' in st_defaults:
            calc['strategy']['aperture_size'] = st_defaults['aperture_size']

    # in most cases sky_key will be the same as ap_key.  however, in the case of IFUNodApPhot, aperture
    # is used but NOT sky_annulus.  this is here to deal with that scenario, but may come up in future
    # scenarios as well.  if sky_annulus is not to be used, then sky_key will be set to None.

    #  NOTE this should get refactored to be more generalized.  use of sky_key for coronagraphy
    #  parameters is hackish...
    if "sky_key" in st_defaults:
        sky_key = st_defaults['sky_key'][calc['strategy']['method']]
    else:
        sky_key = ap_key

    # make sure sky_key is set and sky_annulus_sizes are actually defined
    params = {
        "sky_annulus_sizes": "sky_annulus",
        "unocculted_xys": "unocculted_xy",
        "contrast_azimuths": "contrast_azimuth",
        "contrast_separations": "contrast_separation"
    }
    for p in params:
        if sky_key is not None and p in st_defaults:
            key = calc['configuration']['instrument'][sky_key]
            calc['strategy'][params[p]] = st_defaults[p][key]

    return calc


def build_wb_scene(wb, scene_id=1):
    """
    Take a workbook object and a scene ID (string or int) and build an engine input API
    compatible scene.

    Parameters
    ----------
    wb: dict
        Pandeia workbook format dictionary
    scene_id: int or str (default: 1)
        Key to look up scene definition in workbook

    Returns
    -------
    scene: list of dicts
        Requested scene in engine input API format
    """
    position_keys = ["x_offset", "y_offset", "orientation"]
    try:
        scene_elements = wb['scenes'][str(scene_id)]['scene_elements']
    except KeyError as e:
        msg = "Invalid scene ID or scenes dict: %s (%s)" % (repr(scene_id), repr(e))
        raise EngineInputError(value=msg)
    sources = wb['sources']
    scene = []
    for i, element in scene_elements.items():
        src_id = element['source_id']
        src = sources[str(src_id)]['characteristics']
        src['position'] = {}
        for k in position_keys:
            src['position'][k] = element[k]
        scene.append(src)
    return scene


def build_wb_calc(wb, calc_id=1):
    """
    Take a workbook object and a calculation ID (string or int) and build an engine input API
    compatible calculation.

    Parameters
    ----------
    wb: dict
        Pandeia workbook format dictionary
    calc_id: int or str (default: 1)
        Key to look up calculation definition in workbook

    Returns
    -------
    calc: dict
        Requested calculation in engine input API format
    """
    try:
        calc_input = wb['calculations'][str(calc_id)]
    except KeyError as e:
        msg = "Invalid calculation ID or calculations dict: %s (%s)" % (repr(calc_id), repr(e))
        raise EngineInputError(value=msg)
    calc = {}
    calc['scene'] = build_wb_scene(wb, calc_input['scene_id'])
    if isinstance(calc_input['background'], dict):
        calc['background'] = calc_input['background']['bg_type']
    else:
        calc['background'] = calc_input['background']
    # hack for now to make sure it's something valid. should add hook to BMG here to do it right.
    if calc['background'] == 'dated':
        calc['background'] = 'medium'
    calc['configuration'] = calc_input['camera_config']
    calc['configuration']['instrument']['instrument'] = calc_input['instrument']
    calc['configuration']['instrument']['mode'] = calc_input['insmode']
    calc['strategy'] = calc_input['strategy_args']
    if 'psf_subtraction_source' in calc_input['strategy_args']:
        calc['strategy']['psf_subtraction_source'] = calc['scene'][calc_input['strategy_args']['psf_subtraction_source'] - 1]
    if 'target_type' in calc['strategy'] and 'target_xy' in calc['strategy']:
        # is target_type is source, look up the source, get its position, and configure target_xy accordingly
        if calc['strategy']['target_type'] == 'source':
            src = int(calc['strategy']['target_source']) - 1
            src_cfg = calc['scene'][src]
            xpos = src_cfg['position']['x_offset']
            ypos = src_cfg['position']['y_offset']
            calc['strategy']['target_xy'] = [xpos, ypos]

    return calc


def wb_to_calcs(wb):
    """
    Take a workbook object or filename and create list of engine input API compatible calculations

    Parameters
    ----------
    wb: dict (engine input API) or str
        Workbook data or string containing a workbook filename

    Returns
    -------
    calc_list: list of dicts (engine input API)
        List of dicts in engine input API compatible format
    """
    if isinstance(wb, str):
        try:
            wb = read_json(wb, raise_except=True)
        except:
            msg = "Invalid filename or format for requested workbook %s" % wb
            raise EngineInputError(value=msg)
    calc_list = []
    for i in sorted(wb['calculations']):
        c = build_wb_calc(wb, i)
        calc_list.append(c)
    return calc_list


def calcs_to_wb(calc_list, name="Generated workbook", note="", desc="", wb_id=1, proposal_id="", proposal_state=""):
    """
    Take a list of calculations and make a workbook out of them.  This routine doesn't do any redundancy checking
    between scenes and sources so each calculation has its own scene with its own distinct sources.

    Parameters
    ----------
    calc_list: list of dicts (engine input API)
        List of dicts in engine input API compatible format

    Returns
    -------
    wb: dict
        Workbook data in pandeia workbook format
    """
    wb = {}
    scenes = {}
    sources = {}
    source_i = 1
    calcs = {}
    calc_i = 1
    for c in calc_list:
        # set up the scene for this calculation
        scenes[str(calc_i)] = {}
        scenes[str(calc_i)]["deleted"] = 0
        scenes[str(calc_i)]["characteristics"] = "{}"
        scenes[str(calc_i)]["name"] = "Scene #%d" % calc_i
        scenes[str(calc_i)]["desc"] = "Scene for calculation #%d" % calc_i
        scenes[str(calc_i)]["scene_elements"] = {}
        for src in c['scene']:
            # go through the scene and add the sources there to the workbook's list of sources
            scenes[str(calc_i)]["scene_elements"][str(source_i)] = src['position']
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['name'] = "%d Calculation %s" % (source_i, calc_i)
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['deleted'] = 0
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['source_id'] = source_i
            scenes[str(calc_i)]["scene_elements"][str(source_i)]['source_arg_blob'] = []
            sources[str(source_i)] = {}
            sources[str(source_i)]['deleted'] = 0
            sources[str(source_i)]['name'] = "%f %s %s source" % (
                src['spectrum']['normalization']['norm_flux'],
                src['spectrum']['normalization']['norm_fluxunit'],
                src['shape']['geometry']
            )
            sources[str(source_i)]['characteristics'] = {}
            sources[str(source_i)]['characteristics']['shape'] = src['shape']
            sources[str(source_i)]['characteristics']['spectrum'] = src['spectrum']
            source_i += 1
        # set up the calculations
        calcs[str(calc_i)] = {}
        calcs[str(calc_i)]["strategy_args"] = c['strategy']
        calcs[str(calc_i)]["name"] = "Calculation #%d" % calc_i
        calcs[str(calc_i)]["deleted"] = 0
        calcs[str(calc_i)]["apt"] = 0
        calcs[str(calc_i)]["strategy"] = c['strategy']['method']
        calcs[str(calc_i)]["client_data"] = {
            "cause": "generated_wb"
        }
        calcs[str(calc_i)]["background"] = {"bg_type": c['background'], "ra": 0.0, "dec": 0.0}
        calcs[str(calc_i)]["scene_id"] = calc_i
        calcs[str(calc_i)]["camera_config"] = c['configuration']
        calcs[str(calc_i)]["instrument"] = calcs[str(calc_i)]["camera_config"]["instrument"].pop("instrument")
        calcs[str(calc_i)]["insmode"] = calcs[str(calc_i)]["camera_config"]["instrument"].pop("mode")
        calc_i += 1
    # now populate the workbook
    wb["test_mode"] = "1"
    wb["created"] = 0
    wb["deleted"] = 0
    wb["id"] = wb_id
    wb["name"] = name
    wb["desc"] = desc
    wb["note"] = note
    wb["proposal_id"] = proposal_id
    wb["proposal_state"] = proposal_state
    wb["calculations"] = calcs
    wb["scenes"] = scenes
    wb["sources"] = sources
    return wb
