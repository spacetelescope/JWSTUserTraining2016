#!/usr/bin/env python
"""
Check for required dependencies for the workshop.

Usage::

  % python check_env.py
"""

from distutils.version import LooseVersion


def check_package(package_name, minimum_version=None):
    errors = False
    try:
        pkg = __import__(package_name)
    except ImportError as err:
        print('Error: Failed import: {0}'.format(err))
        errors = True
    else:
        if minimum_version is not None:
            if package_name == 'xlwt':
                installed_version = pkg.__VERSION__
            else:
                installed_version = pkg.__version__
            if (LooseVersion(installed_version) <
                    LooseVersion(str(minimum_version))):
                print('Error: {0} version {1} or later is required, you '
                      'have version {2}'.format(package_name, minimum_version,
                                                installed_version))
                errors = True
    return errors

pkgs = {'IPython': '5.1',
        'notebook': '4.2.3',
        'numpy': '1.6',
        'scipy': '0.15',
        'matplotlib': '1.3',
        'astropy': '1.2.1',
        'photutils': '0.2.2',
        'skimage': '0.12.3',
        'pandas': '0.18.1',
        'glue': '0.8.2',
        'imexam': '0.6',
        'ginga': '2.5.20161005204600'
        }

errors = []
for package_name, min_version in pkgs.items():
    errors.append(check_package(package_name, minimum_version=min_version))
if any(errors):
    print('\nThere are errors that you must resolve before running the '
          'tutorials.')
else:
    print('\nYour Python environment is good to go!')
