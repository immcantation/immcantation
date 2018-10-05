#!/usr/bin/env python3
"""
Build versioning helper script
"""

# Imports
import sys
import yaml
from argparse import ArgumentParser

# Defaults
default_build_file='/Build.yaml'


def readBuild(build_file=default_build_file):
    """
    Read a YAML build file

    Arguments:
      build_file : YAML file containing build information.

    Returns:
      dict : Build file field/value pairs.
    """
    try:
        with open(build_file, 'r') as handle:
            return yaml.load(handle)
    except:
        return None


def writeBuild(field, value, build_file=default_build_file):
    """
    Write build field

    Arguments:
      field : field name.
      value : value for the field.
      build_file : YAML file containing build information.

    Returns:
      True : adds the field/value pair to the file, updating if already present
    """
    # Read existing file
    build = readBuild(build_file)

    # Define default ordering of existing fields
    if build is not None:
        build[field] = value
        order = ['date',
                 'immcantation',
                 'presto',
                 'changeo',
                 'alakazam',
                 'shazam',
                 'tigger',
                 'rdi',
                 'prestor']
        order = [x for x in order if x in build]
    else:
        order = []

    with open(build_file, 'w') as handle:
        for k in order:
            yaml.dump({k: build[k]}, handle, default_flow_style=False)
        if field not in order:
            yaml.dump({field: value}, handle, default_flow_style=False)

    return True


def reportBuild(build_file=default_build_file):
    """
    Report build information

    Arguments:
      build_file : YAML file containing build information.

    Returns:
      str : report.
    """
    # Read build
    build = readBuild(build_file)

    # Set default ordering
    order = ['date',
             'immcantation',
             'presto',
             'changeo',
             'alakazam',
             'shazam',
             'tigger',
             'rdi',
             'prestor']
    order = [x for x in order if x in build]

    report = ['%s: %s' % (k, build[k]) for k in order] + \
             ['%s: %s' % (k, build[k]) for k in build if k not in order]

    print('\n'.join(report))

    return(report)


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
    parser_write = subparsers.add_parser('write',
                                          help='Write fields to a build file.',
                                          description='Write fields to a build file.')
    parser_write.add_argument('-n', action='store', dest='field', type=str, required=True,
                               help='Field name.')
    parser_write.add_argument('-v', action='store', dest='value', type=str, required=True,
                               help='Value for the field.')
    parser_write.add_argument('-f', action='store', dest='build_file', type=str,
                               default=default_build_file, help='YAML build file.')
    parser_write.set_defaults(main=writeBuild)

    # Inspect installed applications
    parser_report = subparsers.add_parser('report',
                                          help='Retrieve build information.',
                                          description='Retrieve build information.')
    parser_report.add_argument('-f', action='store', dest='build_file', type=str,
                               default=default_build_file, help='YAML build file.')
    parser_report.set_defaults(main=reportBuild)

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