#!/usr/bin/env python3
"""
Build pipeline docs from usage statements
"""

# Imports
import os
import re
import subprocess
import yaml
from collections import OrderedDict

# Set YAML loader to OrderedDict
def dict_representer(dumper, data):  return dumper.represent_dict(data.iteritems())
def dict_constructor(loader, node):  return OrderedDict(loader.construct_pairs(node))
yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)

# Load pipeline list
with open(os.path.abspath('../docker/immcantation-release/Pipeline.yaml'), 'r') as f:
    pipelines = yaml.load(f, Loader=yaml.FullLoader)

# Capture usage statements
usage = OrderedDict()
for k, v in pipelines.items():
    x = subprocess.check_output([os.path.abspath('../pipelines/%s' % v[0]), '-h'], text=True)
    usage[k] = re.sub(r'Usage: .*%s' % v[0], r'Usage: %s' % k, x)
    if v[0].endswith('.R'):
        usage[k] = re.sub(r'\n\n', r'\n', usage[k])
        usage[k] = re.sub(r'\nOptions:\n', r'', usage[k])

# Write usage rst file
# Mark for start-after and end-before rst directives
out_path = '_include'
if not os.path.exists(out_path):  os.mkdir(out_path)
with open(os.path.abspath('_include/usage.rst'), 'w') as f:
    for k, v in usage.items():
        f.write('.. Start %s\n\n' % k)
        f.write(v)
        f.write('\n.. End %s\n\n' % k)
