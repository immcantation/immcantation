#!/usr/bin/env python3
"""
Versioning helper script
"""
# Imports
import hglib
import re
import sys
import yaml
from argparse import ArgumentParser
from subprocess import check_output, STDOUT

# Defaults
default_version_file='/Version.yaml'
default_package='immcantation'


class Version():
    """
    Version set class
    """
    def __init__(self, versions):
        self.version = versions['immcantation']['version']
        self.date = versions['immcantation']['date']
        self.packages = {'immcantation': versions['immcantation']['version']}
        self.packages.update(versions['release'])
        self.packages.update(versions['software'])

    def package(self, x):
        return self.packages[x]


def readVersions(version_file=default_version_file):
    """
    Read a YAML version file

    Arguments:
      version_file : YAML file containing version information.

    Returns:
      Version : Version object.
    """
    with open(version_file, 'r') as handle:
        return Version(yaml.load(handle))


def inspectVersions(version_file=default_version_file):
    """
    Determine installed package versions

    Arguments:
      version_file : YAML file containing version information.

    Returns:
      dict : Version strings.
    """
    # Command line version inspection
    #
    # PRESTO=$(python3 -c "import presto; print('%s-%s' % (presto.__version__, presto.__date__))")
    # CHANGEO=$(python3 -c "import changeo; print('%s-%s' % (changeo.__version__, changeo.__date__))")
    # ALAKAZAM=$(Rscript -e "cat(packageDescription('alakazam', fields='Version'))")
    # IGBLAST=$(igblastn -version  | grep 'Package' |sed s/'Package: '//)
    # MUSCLE=$(muscle -version)
    # VSEARCH=$(vsearch --version 2> >(head -n1 | cut -d',' -f1))
    # BLAST=$(blastn -version  | grep 'Package' |sed s/'Package: '//)

    # Only available via the version file
    v = readVersions(version_file=version_file)
    versions = {'immcantation': '%s-%s' % (v.version, v.date)}

    # Python packges
    import presto, changeo, airr
    versions['presto'] = '%s-%s' % (presto.__version__, presto.__date__)
    versions['changeo'] = '%s-%s' % (changeo.__version__, changeo.__date__)
    versions['airr-py'] = '%s' % airr.__version__

    # R packages
    alakazam = check_output('Rscript -e \"cat(packageDescription(\'alakazam\', fields=\'Version\'))\"',
                            stderr=STDOUT, shell=True)
    shazam = check_output('Rscript -e \"cat(packageDescription(\'shazam\', fields=\'Version\'))\"',
                          stderr=STDOUT, shell=True)
    tigger = check_output('Rscript -e \"cat(packageDescription(\'tigger\', fields=\'Version\'))\"',
                          stderr=STDOUT, shell=True)
    rdi = check_output('Rscript -e \"cat(packageDescription(\'rdi\', fields=\'Version\'))\"',
                        stderr=STDOUT, shell=True)
    prestor = check_output('Rscript -e \"cat(packageDescription(\'prestor\', fields=\'Version\'))\"',
                           stderr=STDOUT, shell=True)
    airr_r = check_output('Rscript -e \"cat(packageDescription(\'airr\', fields=\'Version\'))\"',
                           stderr=STDOUT, shell=True)

    versions['alakazam'] = alakazam.decode('utf-8')
    versions['shazam'] = shazam.decode('utf-8')
    versions['tigger'] = tigger.decode('utf-8')
    versions['rdi'] = rdi.decode('utf-8')
    versions['prestor'] = prestor.decode('utf-8')
    versions['airr-r'] = airr_r.decode('utf-8')

    # External applications
    muscle = check_output('muscle -version', stderr=STDOUT, shell=True)
    vsearch = check_output('vsearch --version', stderr=STDOUT, shell=True)
    blast = check_output('blastn -version', stderr=STDOUT, shell=True)
    igblast = check_output('igblastn -version', stderr=STDOUT, shell=True)
    cdhit = check_output('cd-hit-est -h; exit 0', stderr=STDOUT, shell=True)
    phylip = check_output('echo "NULL" | drawtree; exit 0', stderr=STDOUT, shell=True)

    versions['muscle'] = muscle.decode('utf-8').split()[1]
    versions['vsearch'] = re.search(r'(v[0-9.]+)',
                                    vsearch.decode('utf-8').split('\n')[0]).group(0)
    versions['blast'] = re.search(r'(?<=blast )([0-9.]+)',
                                  blast.decode('utf-8').split('\n')[1]).group(0)
    versions['igblast'] = re.search(r'(?<=igblast )([0-9.]+)',
                                    igblast.decode('utf-8').split('\n')[1]).group(0)
    versions['phylip'] =  re.search(r'(?<=PHYLIP version )([0-9.]+)',
                                    phylip.decode('utf-8').split('\n')[0]).group(0)
    versions['cd-hit'] = re.search(r'(?<=CD-HIT version )([0-9.]+)',
                                   cdhit.decode('utf-8').split('\n')[0]).group(0)

    return versions


def updateChangeset(package, repo, version_file):
    """
    Print version for package

    Arguments:
      package : name of the package to return version information for.
      repo : Path to mercurial repository.
      version_file : YAML file containing version information.

    Returns:
      str : changeset updated to.
    """
    # Get version and changeset
    version = getVersion(package, version_file=version_file)
    changeset = getChangeset(version, repo=repo)

    # Update repo
    if changeset is not None:
        client = hglib.open(repo)
        client.update(changeset)

    return(changeset)


def getChangeset(version, repo):
    """
    Print version for package

    Arguments:
      version : Version string to search in tags for.
      repo : Path to mercurial repository.

    Returns:
      str : changeset.
    """
    if version is None:
        print(None)
        return None

    # Build regex
    v = re.compile(r'(^|[\svV])' + version + r'([\s-]|$)')

    # Open repo and retrieve tags
    client = hglib.open(repo)
    tags = client.tags()

    # Check for version number in tags
    changeset = None
    for x in tags:
        if v.search(x[0].decode('utf-8')):
            changeset = '%i:%s' % (x[1], x[2].decode('utf-8'))
            break
    print(changeset)

    return changeset


def getVersion(package=default_package, version_file=default_version_file):
    """
    Print version for package

    Arguments:
      package : name of the package to return version information for.
      version_file : YAML file containing version information.

    Returns:
      str : version.
    """
    v = readVersions(version_file)
    p = v.package(package)

    print(p)
    return(p)


def reportVersions(version_file=default_version_file):
    """
    Report all versions

    Arguments:
      version_file : YAML file containing version information.

    Returns:
      str : version.
    """
    versions = inspectVersions(version_file=version_file)

    release = ['presto',
               'changeo',
               'alakazam',
               'shazam',
               'tigger',
               'rdi',
               'prestor']
    software = ['muscle',
                'vsearch',
                'cd-hit',
                'blast',
                'igblast',
                'phylip',
                'airr-py',
                'airr-r']

    report = ['Immcantation: %s' % versions['immcantation']] + \
             [''] + \
             ['  %s: %s' % (x, versions[x]) for x in release] + \
             [''] + \
             ['  %s: %s' % (x, versions[x]) for x in software]

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

    # Get build file version
    parser_get = subparsers.add_parser('get',
                                       help='Retrieve version number from version file.',
                                       description='Retrieve version number from version file.')
    parser_get.add_argument('-n', action='store', dest='package', type=str, required=True,
                            help='Package name.')
    parser_get.add_argument('-f', action='store', dest='version_file', type=str,
                            default=default_version_file, help='YAML version file.')
    parser_get.set_defaults(main=getVersion)

    # Get mercurial changeset
    parser_cset = subparsers.add_parser('changeset',
                                        help='Retrieve changeset for a tagged version from a mercurial repository.',
                                        description='Retrieve changeset for a tagged version from a mercurial repository.')
    parser_cset.add_argument('-v', action='store', dest='version', type=str, required=True,
                             help='Version number.')
    parser_cset.add_argument('-r', action='store', dest='repo', type=str, required=True,
                             help='Path to the mercurial repository.')
    parser_cset.set_defaults(main=getChangeset)

    # Update mercurial changeset
    parser_update = subparsers.add_parser('update',
                                          help='Update a mercurial repository to tagged version number from version file.',
                                          description='Update a mercurial repository to tagged version number from version file.')
    parser_update.add_argument('-n', action='store', dest='package', type=str, required=True,
                               help='Package name.')
    parser_update.add_argument('-r', action='store', dest='repo', type=str, required=True,
                               help='Path to the mercurial repository.')
    parser_update.add_argument('-f', action='store', dest='version_file', type=str,
                               default=default_version_file, help='YAML version file.')
    parser_update.set_defaults(main=updateChangeset)

    # Inspect installed applications
    parser_report = subparsers.add_parser('report',
                                          help='Retrieve version information from installed packages.',
                                          description='Retrieve version information from installed packages.')
    parser_report.add_argument('-f', action='store', dest='version_file', type=str,
                               default=default_version_file, help='YAML version file.')
    parser_report.set_defaults(main=reportVersions)

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