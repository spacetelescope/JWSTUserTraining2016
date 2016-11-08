
# User Training in JWST Data Analysis II

The second ["User Training in JWST Data
Analysis"](https://jwst.stsci.edu/events/events-area/stsci-events-listing-container/user-training-in-jwst-data-analysis-ii)
workshop will be held on November 8-11, 2016 at the Space Telescope
Science Institute (STScI).  The purpose of this three-day meeting is
to provide training in the the open-source data analysis tools being
developed at STScI for use with the James Webb Space Telescope (JWST),
as well as for many other optical/IR observatories.

The meeting will include an optional training day (Day Zero) on
November 8 to provide novice participants an introduction to Python
and Astropy.

## Installation/Setup

This will be an interactive workshop so be sure to come with a laptop
prepared to try out some of the tools that will be discussed and
demoed.  Before arriving at the workshop next week, we ask that you
install the the Anaconda distribution for Python 3.5, which we have
packaged along with some additional software.  Downloads for Mac and
Linux can be found at:

* http://ssb.stsci.edu/conda/installers/AstroConda-1.0.2-Linux-x86_64.sh
* http://ssb.stsci.edu/conda/installers/AstroConda-1.0.2-MacOSX-x86_64.sh

These include both Anaconda and the AstroConda software repository,
which contains additional tools that will be shown at the workshop.
If you have trouble installing using the above files (or are using
Windows), you will need to download anaconda separately
(https://www.continuum.io/downloads), and then install Astroconda on
top following the instructions here:
http://astroconda.readthedocs.io/en/latest/ .

If you already have Anaconda installed on your machine, you can create
a special environment for this workshop which contains all the software
you will need using this environment file:

```shell
% conda env create -n jwst-workshop --file environment.yml
```

The command above will create an environment called "jwst-workshop",
but you can change that to any other desirable name.


You can run the ``check_env.py`` script to perform a basic check of your
Python environment and some of the required dependencies::

```shell
% python check_env.py
```

If you have issues getting set up, you can also run the notebooks on mybinder.org:

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/spacetelescope/jwstusertraining2016)
