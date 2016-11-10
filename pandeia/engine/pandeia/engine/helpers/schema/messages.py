#
# Full schema for basic all-encompassing workbook checking.
#
# Note the use of the additionalProperties keyword. Because our
# collections (e.g. calculations, scenes, sources) are dicts and
# not arrays, with the key being the index number, to do this schema
# the simplest way, we would have multiple lines in each collection,
# one for each index (1-N), which would be silly.
#
# For the use of "oneOf", this page was very helpful:
#    http://json-schema.org/example2.html
# see also:
#    http://spacetelescope.github.io/understanding-json-schema
#
# This module takes advantage of Python's dictionaries (and sub-dicts)
# to keep things modular, and to create (below) schemas that allow
# validation of a whole workbook, or just of parts (like a source).
#
# Keep in synch with:
#     pandeia/engine/doc/engine_input_api.rst
#

# Globals
SCHEMA_URL = 'http://json-schema.org/schema#'

# an empty dict ({}) is valid for this schema but not an empty string('')
EMPTY_SCHEMA = {
    'type': 'object',
    '$schema': SCHEMA_URL,
    'additionalProperties': False
}

# This regex matches a string representation of an integer.
INT_REGEX = '^-?[0-9]+$'

#-------------------------------------------------------------------------------
# DEFINITIONS OF PARTS - these are not to be used directly during validation,
# but instead incorporated into regular schemas.
#
# These are declared in order of most atomic to most complex
#-------------------------------------------------------------------------------

a_geom_point = {
    'type': 'object',
    'required': ['geometry'],
    'additionalProperties': False,
    'properties': {
        'geometry': {'type': 'string', 'enum': ['point']},
    },
}

a_geom_2d = {
    'type': 'object',
    'required': ['geometry','major','minor'],
    'additionalProperties': False,
    'properties': {
        'geometry': {'type': 'string', 'enum': ['flat','gaussian2d']},
        'major':    {'type': 'number'},
        'minor':    {'type': 'number'},
    },
}

a_geom_sersic = {
    'type': 'object',
    'required': ['geometry','major','minor','sersic_index'],
    'additionalProperties': False,
    'properties': {
        'geometry': {'type': 'string', 'enum': ['sersic']},
        'major':        {'type': 'number'},
        'minor':        {'type': 'number'},
        'sersic_index': {'type': 'number'},
    },
}

a_normalization_none = {
    'type': 'object',
    'required': ['type'],
    'additionalProperties': False,
    'properties': {
        'type': {'type': 'string', 'enum': ['none']},
    },
}

a_normalization_at_lambda = {
    'type': 'object',
    'required': ['norm_flux','norm_fluxunit','norm_wave','norm_waveunit','type'],
    'additionalProperties': False,
    'properties': {
        'norm_flux':     {'type': 'number'},
        'norm_fluxunit': {'type': 'string', 'enum': ['flam','fnu','abmag','vegamag','mjy','ujy','njy','jy']},
        'norm_wave':     {'type': 'number'},
        'norm_waveunit': {'type': 'string', 'enum': ['microns']},
        'type':          {'type': 'string', 'enum': ['at_lambda']},
    },
}

a_normalization_instrument = {
    'type': 'object',
    'required': ['bandpass','norm_flux','norm_fluxunit','type'],
    'additionalProperties': False,
    'properties': {
        'bandpass':      {'type': 'string'},
        'norm_flux':     {'type': 'number'},
        'norm_fluxunit': {'type': 'string', 'enum': ['flam','fnu','abmag','vegamag','mjy','ujy','njy','jy']},
        'type':          {'type': 'string', 'enum': ['hst','jwst','photsys']},
    },
}

an_sed_none = {
    'type': 'object',
    'required': ['sed_type'],
    'additionalProperties': False,
    'properties': {
        'sed_type': {'type': 'string', 'enum': ['no_continuum']},
    },
}

an_sed_flat = {
    'type': 'object',
    'required': ['sed_type','unit'],
    'additionalProperties': False,
    'properties': {
        'sed_type': {'type': 'string', 'enum': ['flat']},
        'unit':     {'type': 'string', 'enum': ['flam','fnu'] },
    },
}

an_sed_blackbody = {
    'type': 'object',
    'required': ['sed_type','temp'],
    'additionalProperties': False,
    'properties': {
        'sed_type': {'type': 'string', 'enum': ['blackbody']},
        'temp':     {'type': 'number'},
    },
}

an_sed_powerlaw = {
    'type': 'object',
    'required': ['index','sed_type','unit'],
    'additionalProperties': False,
    'properties': {
        'index':    {'type': 'number'},
        'sed_type': {'type': 'string', 'enum': ['powerlaw']},
        'unit':     {'type': 'string', 'enum': ['flam','fnu'] },
    },
}

an_sed_keyed = {
    'type': 'object',
    'required': ['key','sed_type'],
    'additionalProperties': False,
    'properties': {
        'key':      {'type': 'string'},
        'sed_type': {'type': 'string', 'enum': ['brown','hst_calspec','phoenix']},
    },
}

an_sed_upload = {
    'type': 'object',
    'required': ['sed_type','spectrum_id'],
    'additionalProperties': False,
    'properties': {
        'sed_type':    {'type': 'string', 'enum': ['input']},
        'spectrum_id': {'type': 'number'},
    },
}

a_background_none = {
    'type': 'object',
    'required': ['bg_type'],
    'additionalProperties': False,
    'properties': {
        'bg_type': {'type': 'string', 'enum': ['none']},
    },
}

a_background_dateless = {
    'type': 'object',
    'required': ['bg_type','ra','dec','ra_dec_str'],
    'additionalProperties': False,
    'properties': {
        'bg_type':    {'type': 'string', 'enum': ['low','medium','high']},
        'ra':         {'type': 'number'},
        'dec':        {'type': 'number'},
        'ra_dec_str': {'type': 'string'},
    },
}

a_background_dated = {
    'type': 'object',
    'required': ['bg_type','ra','dec','date'],
    'additionalProperties': False,
    'properties': {
        'bg_type':    {'type': 'string', 'enum': ['low','medium','high','dated']},
        'ra':         {'type': 'number'},
        'dec':        {'type': 'number'},
        'date':       {'type': 'number'},
        'ra_dec_str': {'type': 'string'}, # optional
        'date_str':   {'type': 'string'}, # optional
    },
}

a_background_choice = {
    'type': 'object',
    'oneOf': [a_background_none, a_background_dateless, a_background_dated],
}

a_helper_strategy_background = {
    'type': 'object',
    'required': ["background_subtraction"],

    'properties': {
        'background_subtraction':   {'type': 'boolean', },
        'sky_annulus':   {'type': 'array', 'items': {'type': 'number'}},
    },

    # Dep from background_subtraction to sky_annulus
    'oneOf': [
        {
            "properties": {
                'background_subtraction': { 'enum': [True] },
            },
            "required": ["sky_annulus"],
        },
        {
            "properties": {
                'background_subtraction': { 'enum': [False] },
            },
        },
    ],

    # Dep from sky_annulus to background_subtraction
    'dependencies': {
        'sky_annulus': {
            'properties': {
                'background_subtraction': { 'enum': [True] },
            },
        },
    },
}

a_strategy_coronagraphy = {
    'allOf': [
        {
            'type': 'object',
            'required': ['annulus_shape','aperture_size','contrast_azimuth','contrast_separation','method','psf_subtraction','psf_subtraction_source','scene_rotation','target_source','target_type','target_xy','units'],
            'properties': {
                'annulus_shape':          {'type': 'string', 'enum': ['circular']},
                'aperture_size':          {'type': 'number'},
                'background_subtraction': {'type': 'boolean', 'enum': [True]},
                'contrast_azimuth':       {'type': 'number'},
                'contrast_separation':    {'type': 'number'},
                'method':                 {'type': 'string', 'enum': ['coronagraphy']},
                'psf_subtraction':        {'type': 'string', 'enum': ['none','realistic','optimal']},
                'psf_subtraction_source': {'type': 'number'},
                'scene_rotation':         {'type': 'number'},
                'target_source':          {'type': 'string', 'pattern': INT_REGEX},
                'target_type':            {'type': 'string', 'enum': ['source']},
                'target_xy':              {'type': 'array', 'items': {'type': 'number'}},
                'ui_strat_xy':            {'type': 'array', 'items': {'type': 'number'}}, # optional
                'units':                  {'type': 'string', 'enum': ['arcsec']},
            },
        },
        a_helper_strategy_background,
    ],
}

a_strategy_ifunodinscene = {
    'type': 'object',
    'required': ['aperture_size','dithers','method','reference_wavelength','target_source','target_type','target_xy','units'],
    'additionalProperties': False,
    'properties': {
        'aperture_size':        {'type': 'number'},
        'dithers': {
            'type': 'array',
            'items': {
                'type': 'object',
                'required': ['x','y'],
                'additionalProperties': False,
                'properties': {
                    'x': {'type': 'number'},
                    'y': {'type': 'number'},
                }
            }
        },
        'method':               {'type': 'string', 'enum': ['ifunodinscene']},
        'reference_wavelength': {'oneOf': [{'type': 'number'}, {'type': 'null'}]},
        'target_source':        {'type': 'string', 'pattern': INT_REGEX},
        'target_type':          {'type': 'string', 'enum': ['coords','source']},
        'target_xy':            {'type': 'array', 'items': {'type': 'number'}},
        'ui_strat_xy':          {'type': 'array', 'items': {'type': 'number'}}, # optional
        'units':                {'type': 'string', 'enum': ['arcsec']},
    },
}

a_strategy_ifunodoffscene = {
    'type': 'object',
    'required': ['aperture_size','method','reference_wavelength','target_source','target_type','target_xy','units'],
    'additionalProperties': False,
    'properties': {
        'aperture_size':        {'type': 'number'},
        'dithers': { # optional for UI, used only by engine
            'type': 'array',
            'items': {
                'type': 'object',
                'required': ['x','y'],
                'additionalProperties': False,
                'properties': {
                    'x': {'type': 'number'},
                    'y': {'type': 'number'},
                }
            }
        },
        'method':               {'type': 'string', 'enum': ['ifunodoffscene']},
        'reference_wavelength': {'oneOf': [{'type': 'number'}, {'type': 'null'}]},
        'target_source':        {'type': 'string', 'pattern': INT_REGEX},
        'target_type':          {'type': 'string', 'enum': ['coords','source']},
        'target_xy':            {'type': 'array', 'items': {'type': 'number'}},
        'ui_strat_xy':          {'type': 'array', 'items': {'type': 'number'}}, # optional
        'units':                {'type': 'string', 'enum': ['arcsec']},
    },
}

a_strategy_imagingapphot = {
    'allOf': [
        {
            'type': 'object',
            'required': ['aperture_size','method','target_source','target_type','target_xy','units'],
            'properties': {
                'aperture_size': {'type': 'number'},
                'method':        {'type': 'string', 'enum': ['imagingapphot']},
                'target_source': {'type': 'string', 'pattern': INT_REGEX},
                'target_type':   {'type': 'string', 'enum': ['coords','source']},
                'target_xy':     {'type': 'array', 'items': {'type': 'number'}},
                'ui_strat_xy':   {'type': 'array', 'items': {'type': 'number'}}, # optional
                'units':         {'type': 'string', 'enum': ['arcsec']},
            },
        },
        a_helper_strategy_background,
    ],
}

a_strategy_msafullapphot = {
    'type': 'object',
    'required': ['background_subtraction','dithers','method','reference_wavelength','shutter_offset','target_source','target_type','units'],
    'additionalProperties': False,
    'properties': {
        'background_subtraction':   {'type': 'boolean', },
        'dithers': {
            'type': 'array',
            'items': {
                'type': 'object',
                'required': ['x','y','on_source'],
                'additionalProperties': False,
                'properties': {
                    'x':         {'type': 'number'},
                    'y':         {'type': 'number'},
                    'on_source': {'type': 'array', 'items': {'type': 'boolean'}},
                },
            },
        },
        'method':               {'type': 'string', 'enum': ['msafullapphot']},
        'reference_wavelength': {'oneOf': [{'type': 'number'}, {'type': 'null'}]},
        'shutter_offset':       {'type': 'array', 'items': {'type': 'number'}},
        'target_source':        {'type': 'string', 'pattern': INT_REGEX},
        'target_type':          {'type': 'string', 'enum': ['source']},
        'ui_strat_xy':          {'type': 'array', 'items': {'type': 'number'}}, # optional
        'units':                {'type': 'string', 'enum': ['arcsec']},
    },
}

a_strategy_soss = {
    'type': 'object',
    'required': ['background_subtraction','method','order','reference_wavelength'], # aperture_size, units
    'additionalProperties': False,
    'properties': {
        'background_subtraction': {'type': 'boolean', 'enum': [False]},
        'method':                 {'type': 'string', 'enum': ['soss']},
        'order':                  {'type': 'number'},
        'reference_wavelength':   {'type': 'number'},
    },
}

a_strategy_specapphot = {
    'allOf': [
        {
            'type': 'object',
            'required': ['aperture_size','method','reference_wavelength','target_source','target_type','target_xy','units'],
            'properties': {
                'aperture_size':        {'type': 'number'},
                'method':               {'type': 'string', 'enum': ['specapphot']},
                'reference_wavelength': {'oneOf': [{'type': 'number'}, {'type': 'null'}]},
                'target_source':        {'type': 'string', 'pattern': INT_REGEX},
                'target_type':          {'type': 'string', 'enum': ['coords','source']},
                'target_xy':            {'type': 'array', 'items': {'type': 'number'}},
                'ui_strat_xy':          {'type': 'array'}, # optional (also, vals may be None or float)
                'units':                {'type': 'string', 'enum': ['arcsec']},
            },
        },
        a_helper_strategy_background,
    ],
}

a_strategy_choice = {
    'type': 'object',
    'oneOf': [a_strategy_coronagraphy, a_strategy_ifunodinscene, a_strategy_ifunodoffscene,
              a_strategy_imagingapphot, a_strategy_msafullapphot, a_strategy_soss, a_strategy_specapphot],
}

a_scene_element = {
    'type': 'object',
    'required': ['deleted','name','orientation','source_arg_blob','source_id','x_offset','y_offset'],
    'additionalProperties': False,
    'properties': {
        'deleted':         {'type': 'number'},
        'name':            {'type': 'string'},
        'orientation':     {'type': 'number'},
        'source_arg_blob': EMPTY_SCHEMA, # unused at current
        'source_id':       {'type': 'number'},
        'x_offset':        {'type': 'number'},
        'y_offset':        {'type': 'number'},
    },
}

a_scene = {
    'type': 'object',
    'required': ['characteristics','deleted','desc','name','scene_elements'],
    'additionalProperties': False,
    'properties': {
        'characteristics': EMPTY_SCHEMA, # unused at current
        'deleted':         {'type': 'number'},
        'desc':            {'type': 'string'},
        'name':            {'type': 'string'},
        'scene_elements': {
            'type': 'object',
            'additionalProperties': a_scene_element,
            'properties': {},
        },
    },
}

a_source_line = {
    'type': 'object',
    'required': ['center','emission_or_absorption','name','profile','strength','width'],
    'additionalProperties': False,
    'properties': {
        'center':                 {'type': 'number'},
        'emission_or_absorption': {'type': 'string', 'enum': ['emission','absorption']},
        'name':                   {'type': 'string'},
        'profile':                {'type': 'string', 'enum': ['gaussian']},
        'strength':               {'type': 'number'},
        'width':                  {'type': 'number'},
    },
}

a_source_characteristics = {
    'type': 'object',
    'required': ['shape','spectrum'],
    'additionalProperties': False,
    'properties': {
        'shape': {
            'type': 'object',
            'oneOf': [a_geom_point, a_geom_2d, a_geom_sersic],
        },
        'spectrum': {
            'type': 'object',
            'required': ['lines','normalization','redshift','sed','extinction'],
            'additionalProperties': False,
            'properties': {
                'lines': {'type': 'array', 'items': a_source_line},
                'normalization': {
                    'type': 'object',
                    'oneOf': [a_normalization_none, a_normalization_at_lambda, a_normalization_instrument],
                },
                'redshift': {'type': 'number'},
                'sed': {
                    'type': 'object',
                    'oneOf': [an_sed_none, an_sed_blackbody, an_sed_flat, an_sed_powerlaw, an_sed_keyed, an_sed_upload],
                },
                'extinction': {
                    'law': {'type': 'string'},
                    'bandpass': {'type': 'string'},
                    'unit': {'type': 'string'},
                    'value': {'type': 'float'}
                }
            }
        }
    }
}

a_source = {
    'type': 'object',
    'required': ['characteristics','deleted','name'],
    'additionalProperties': False,
    'properties': {
        'characteristics': a_source_characteristics,
        'deleted': {'type': 'number'},
        'name':    {'type': 'string'},
    },
}

a_camera_config = {
    'type': 'object',
    'required': ['detector','instrument'],
    'additionalProperties': False,
    'properties': {
        'detector': {
            'type': 'object',
            'required': ['nexp','ngroup','nint','readmode'],
            'additionalProperties': False,
            'properties': {
                'nexp':     {'type': 'number'},
                'ngroup':   {'type': 'number'},
                'nint':     {'type': 'number'},
                'readmode': {'type': 'string'},
                'subarray': {'type': 'string'}, # optional
            },
        },
        'instrument': {
            'type': 'object',
            'required': ['aperture','disperser','filter'],
            'additionalProperties': False,
            'properties': {
                'aperture':         {'type': 'string'},
                'disperser':        {'oneOf': [{'type': 'string'}, {'type': 'null'}]},
                'filter':           {'oneOf': [{'type': 'string'}, {'type': 'null'}]},
                'shutter_location': {'type': 'string'}, # optional
                'slitlet_shape': { # optional
                    'type': 'array',
                    'items': {
                        'type': 'array',
                        'items': {'type': 'number'}
                    }
                }
            }
        }
    }
}

a_client_data = {
    'type': 'object' # used only by client - maybe define more explicitly later
}

a_calculation = {
    'type': 'object',
    'required': ['apt','background','camera_config','deleted','insmode','instrument','name','scene_id','strategy','strategy_args'],
    'additionalProperties': False,
    'properties': {
        'apt':             {'type': 'string'},
        'background':      a_background_choice,
        'camera_config':   a_camera_config,
        'client_data':     a_client_data, # optional
        'deleted':         {'type': 'number'},
        'insmode':         {'type': 'string'},
        'instrument':      {'type': 'string', 'enum': ['miri','nircam','niriss','nirspec']},
        'name':            {'type': 'string'},
        'scene_id':        {'type': 'number'},
        'strategy':        {'type': 'string'},
        'strategy_args':   a_strategy_choice,
    }
}


#-------------------------------------------------------------------------------
# Complete, actionable schemas to be used for validation
#-------------------------------------------------------------------------------

WORKBOOK_SCHEMA = {
    '$schema': SCHEMA_URL,
    'type': 'object',
    'required': ['calculations','created','deleted','desc','hash','id','name','note',
                 'proposal_id','proposal_state','scenes','sources','test_mode'],
    'additionalProperties': False,
    'properties': {
        'created':        {'type': 'number'},
        'deleted':        {'type': 'number'},
        'desc':           {'type': 'string'},
        'hash':           {'type': 'string'},
        'id':             {'type': 'number'},
        'name':           {'type': 'string'},
        'note':           {'type': 'string'},
        'proposal_id':    {'type': 'string'},
        'proposal_state': {'type': 'string'},
        'test_mode':      {'type': 'string'},
        'calculations': {
            'type': 'object',
            'additionalProperties': {'$ref': '#/definitions/a_calculation'},
            'properties': {},
        },
        'scenes': {
            'type': 'object',
            'additionalProperties': a_scene,
            'properties': {},
        },
        'sources': {
            'type': 'object',
            'additionalProperties': a_source,
            'properties': {},
        },
    },

    # definitions: anything listed here can be used in the $ref syntax
    'definitions': {

        # definitions from dicts defined above, in order of increasing complexity
        'a_geom_point':                 a_geom_point,
        'a_geom_2d':                    a_geom_2d,
        'a_geom_sersic':                a_geom_sersic,
        'a_normalization_none':         a_normalization_none,
        'a_normalization_at_lambda':    a_normalization_at_lambda,
        'a_normalization_instrument':   a_normalization_instrument,
        'an_sed_none':                  an_sed_none,
        'an_sed_flat':                  an_sed_flat,
        'an_sed_blackbody':             an_sed_blackbody,
        'an_sed_powerlaw':              an_sed_powerlaw,
        'an_sed_keyed':                 an_sed_keyed,
        'a_background_none':            a_background_none,
        'a_background_dateless':        a_background_dateless,
        'a_background_dated':           a_background_dated,
        'a_background_choice':          a_background_choice,

        # strategies
        'a_helper_strategy_background': a_helper_strategy_background,
        'a_strategy_coronagraphy':      a_strategy_coronagraphy,
        'a_strategy_ifunodinscene':     a_strategy_ifunodinscene,
        'a_strategy_ifunodoffscene':    a_strategy_ifunodoffscene,
        'a_strategy_imagingapphot':     a_strategy_imagingapphot,
        'a_strategy_msafullapphot':     a_strategy_msafullapphot,
        'a_strategy_soss':              a_strategy_soss,
        'a_strategy_specapphot':        a_strategy_specapphot,
        'a_strategy_choice':            a_strategy_choice,

        # scenes/sources
        'a_scene_element':              a_scene_element,
        'a_scene':                      a_scene,
        'a_source':                     a_source,

        # calculations
        'a_calculation':                a_calculation,
    },
}

# These are the constant schema objects intended for general use.
# For the following objects, use our dict but add schema entry (shallow copy is good enough)

SOURCE_CHAR_SCHEMA = a_source_characteristics.copy()
SOURCE_CHAR_SCHEMA['$schema'] = SCHEMA_URL

CALC_BG_SCHEMA = a_background_choice.copy()
CALC_BG_SCHEMA['$schema'] = SCHEMA_URL

CALC_STRATEGY_SCHEMA = a_strategy_choice.copy()
CALC_STRATEGY_SCHEMA['$schema'] = SCHEMA_URL

CALC_CAMERA_SCHEMA = a_camera_config.copy()
CALC_CAMERA_SCHEMA['$schema'] = SCHEMA_URL

CALC_CLIENT_DATA_SCHEMA = a_client_data.copy()
CALC_CLIENT_DATA_SCHEMA['$schema'] = SCHEMA_URL

