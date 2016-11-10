# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Module that implements a single calculation API for the ETC3D engine
"""

from __future__ import division, absolute_import

from .etc3D import calculate_exposure_time, calculate_sn


def perform_calculation(calc_input, reverse=False, dict_report=True, webapp=True):
    """
    Function to perform a single ETC calculation

    Parameters
    ----------
    input: dict
        Dictionary containing the information required to perform the calculation.
    reverse: boolean (default: False)
        Specify whether calculation is 'reverse' where a desired signal/noise is specified
        and the calculation determines an optimal ExposureSpecification to achieve it.
    dict_report: Boolean (default: True)
        If True, return a dict in engine output API format. Otherwise return
        a report.Report instance.
    webapp: Boolean (default: True)
        Toggle strict API checking for webapp
    """

    if True:   # you can change it to False
        if 'fake_exception' in calc_input:
            perform_fake_exceptions(calc_input)

    # run the calculation....
    if reverse:
        report = calculate_exposure_time(calc_input, webapp=webapp)
    else:
        report = calculate_sn(calc_input, webapp=webapp)

    if dict_report:
        return report.as_dict()
    else:
        return report


def perform_fake_exceptions(calc_input):
    i = calc_input['fake_exception']
    if 'pyetc' in i:
        from . import custom_exceptions
        raise custom_exceptions.PyetcException("fake pyetc exception for testing")
    if 'exception' in i:
        raise Exception("fake abnormal exception for testing")
