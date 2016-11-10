import urllib2
import json

from pandeia.engine.helpers.accessor_globals import bmg_url
from pandeia.engine.custom_exceptions import BMGError


def get_in_field_bg(ra, dec, date, ra_dec_str, date_str, components=11):

    '''
    component is a bitmap of modes:

        zodi = 01
        ism = 02
        cib = 08

        combined = 11
    '''
    data = {
        'query_type': 'in_field',
        'ra_dec': [ra, dec],
        'date': date,
        'components': components,
        'request_id': '1'
    }

    data_json = json.dumps(data)
    request = '%s/bmg/bmgws/?bmg=%s' % (bmg_url, data_json)
    request = request.replace(' ', '')

    print()
    print('BMG REQUEST:')
    print(request)
    print()

    try:
        response = urllib2.urlopen(request)
    except urllib2.HTTPError as e:
        # handle HTTP errors
        http_status = str(e.code)
        raise BMGError(http_status, ra_dec_str, date_str)

    response = response.read()
    response = json.loads(response)

    status = response['status']
    if status != 0:
        raise BMGError('Non-zero status code returned from BMG for request: '+request, ra_dec_str, date_str)
    else:
        [wave, bmg_bg] = response['spectrum']
        return wave, bmg_bg


def get_zodi(ra, dec, date, ra_dec_str, date_str):
    """
    Wrap get_in_field_bg() to only return the zodiacal background component by using "components=1"
    """
    return get_in_field_bg(ra, dec, date, ra_dec_str, date_str, components=1)


def get_ism(ra, dec, date, ra_dec_str, date_str):
    """
    Wrap get_in_field_bg() to only return the ISM (aka cirrus) background component by using "components=2"
    """
    return get_in_field_bg(ra, dec, date, ra_dec_str, date_str, components=2)


def get_cib(ra, dec, date, ra_dec_str, date_str):
    """
    Wrap get_in_field_bg() to only return the cosmic infrared background component by using "components=8"
    """
    return get_in_field_bg(ra, dec, date, ra_dec_str, date_str, components=8)
