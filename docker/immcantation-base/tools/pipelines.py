#!/usr/bin/env python3
"""
Build pipeline setup helper script
"""

# Imports
import os
import shutil
import sys
import yaml
from argparse import ArgumentParser
from collections import OrderedDict

# Defaults
default_pipelines='/Pipeline.yaml'
default_source='/tmp/immcantation/pipelines'
default_target='/usr/local/bin'

# Set YAML loader to OrderedDict
def dict_representer(dumper, data):  return dumper.represent_dict(data.iteritems())
def dict_constructor(loader, node):  return OrderedDict(loader.construct_pairs(node))
yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)


def readPipelines(pipelines=default_pipelines):
    """
    Read a YAML pipeline file

    Arguments:
      pipelines (str): YAML file containing pipeline information.

    Returns:
      dict: build file field/value pairs.
    """
    with open(pipelines, 'r') as handle:
        return yaml.load(handle, Loader=yaml.FullLoader)


def copyPipelines(pipelines=default_pipelines,
                  source=default_source,
                  target=default_target):
    """
    Copy pipeline scripts

    Arguments:
      pipelines (str): YAML file containing pipeline information.
      source (str): source directory.
      target (str): target directory.
    """
    # Get pipelines
    scripts = readPipelines(pipelines)

    # Define default ordering of existing fields
    for k, v in scripts.items():
        shutil.copy(os.path.join(source, v[0]), os.path.join(target, k))

    return True


def reportPipelines(pipelines=default_pipelines):
    """
    Report pipeline scripts

    Arguments:
      pipelines (str): YAML file containing pipeline information.
    """
    # Get pipelines
    scripts = readPipelines(pipelines)

    # Define default ordering of existing fields
    for k, v in scripts.items():
        print('%s: %s' % (k, v[1]))

    return True


def getArgParser():
    """
    Defines the ArgumentParser

    Returns:
     argparse.ArgumentParser : argument parser
    """
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(title='subcommands', metavar='', help='Task')
    subparsers.required = True

    # Write build information
    parser_write = subparsers.add_parser('copy',
                                          help='Copy pipeline scripts.',
                                          description='Copy pipeline scripts.')
    parser_write.add_argument('-s', action='store', dest='source', type=str,
                              default=default_source, help='Source directory.')
    parser_write.add_argument('-t', action='store', dest='target', type=str,
                              default=default_target, help='Target directory.')
    parser_write.add_argument('-f', action='store', dest='pipelines', type=str,
                              default=default_pipelines, help='YAML pipeline file.')
    parser_write.set_defaults(main=copyPipelines)

    # Inspect installed applications
    parser_report = subparsers.add_parser('report',
                                          help='Retrieve pipeline information.',
                                          description='Retrieve pipeline information.')
    parser_report.add_argument('-f', action='store', dest='pipelines', type=str,
                              default=default_pipelines, help='YAML pipeline file.')
    parser_report.set_defaults(main=reportPipelines)

    return(parser)


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    parser = getArgParser()
    args = parser.parse_args()

    main = args.main
    args_dict = args.__dict__
    del args_dict['main']

    check = main(**args_dict)
    if check is None:  sys.exit(1)