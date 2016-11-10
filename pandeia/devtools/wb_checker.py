#!/usr/bin/env python
#
# Basic all-encompassing workbook checker.  For now it
# starts with only verifying the wb format is correct.
#

import json, jsonschema, os, sys
from pandeia.engine.helpers.schema.messages import WORKBOOK_SCHEMA


def read_full_file(fname):
    """ convenience function """
    assert os.path.exists(fname), 'read_full_file expected, but did not find file named: '+fname
    f = open(fname, 'r')
    buf = f.read()
    f.close()
    return buf

def json_fname_to_dict(fname):
    """ convenience function """
    buf = read_full_file(fname)
    return json.loads(read_full_file(fname))

def wb_has_valid_format(wbobj):
    retval = ''
    try:
        jsonschema.validate(wbobj, WORKBOOK_SCHEMA)
    except Exception as e:
        retval = str(e)
    return retval

def check_wb(wbobj):
    """ Run all known checks on a wb. Returns empty string
        if all is well, otherwise returns all caught errors.
    """
    # for now we only check format
    return wb_has_valid_format(wbobj)


#
# main routine
#
if __name__=='__main__': # in case something else imports this file

    # run checks on all files in arg list
    exitval = 0
    for arg_fname in sys.argv[1:]:
        fname = os.path.abspath(arg_fname)
        print('Checking wb: '+fname)
        obj = json_fname_to_dict(fname)
        res = check_wb(obj)
        if res != '':
            print(res)
            exitval += 1
    sys.exit(exitval)
