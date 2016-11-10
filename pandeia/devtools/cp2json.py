#!/usr/bin/env python

import json
import sys
from collections import defaultdict
import ConfigParser as cp

def nested_defaultdict():
    return defaultdict(nested_defaultdict)

pars = cp.SafeConfigParser()
pars.read(sys.argv[1])

config = nested_defaultdict()

for s in pars.sections():
    key = s.lower()
    if ':' in key:
        keys = key.split(':')
        dict_string = 'config'
        for k in keys:
            dict_string += "['%s']" % k
        sec_data = dict(pars.items(s))
        exec('%s = sec_data' % dict_string)
    else:
        config[key].update(dict(pars.items(s)))

#print repr(config)

with open("output.json", 'w') as f:
    json.dump(config, f, indent=4)

# for k in config.keys():
#     file = "%s.cfg" % k
#     with open(file, 'w') as f:
#         json.dump(config[k], f, indent=4)
