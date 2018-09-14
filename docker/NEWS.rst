Docker Container Release Notes
========================================================================

Version 2.1.0:  September 12, 2018
------------------------------------------------------------------------

Version Updates:

+ alakazam 0.2.11


Version 2.0.0:  September 8, 2018
------------------------------------------------------------------------

Version Updates:

+ pRESTO 0.5.9
+ Change-O 0.4.2
+ airr 1.2.0

Image Changes:

+ Added tbl2asn.

Pipeline Changes:

+ Changed behavior of subsampling argument to ``shazam-threshold``
  to subsample distances after nearest-neighbor distance calculation
  rather than rows before distance calculation.


Version 1.10.2:  July 3, 2018
------------------------------------------------------------------------

Pipeline Changes:

+ Added data set subsampling to ``shazam-threshold`` with a default
  value of 15000 records.
+ Added ``-f`` argument to ``changeo-igblast`` to allow optional
  filtering of non-productive/non-functional sequences.
+ Added ``-a`` argument to ``changeo-clone`` to allow retention of
  non-productive/non-functionals sequences during cloning.
+ Added ``-v`` argument to ``tigger-genotype`` to allow specification of
  the V genotyped column name.


Version 1.10.1:  July 1, 2018
------------------------------------------------------------------------

Pipeline Changes:

+ Fixed a bug wherein ``changeo-igblast`` and ``changeo-clone`` were
  not working with an unspecified output directory (``-o`` argument).
+ Updated CPU core detection in ``tigger-genotype`` and
  ``shazam-threshold`` for compatability with new R package versions.

Accessory Script Changes:

+ Fixed ``fetch_imgtdb.sh`` creating empty mouse IGKC and IGLC files.

Image Changes:

+ Changed default CRAN mirror setting.


Version 1.10.0:  May 23, 2018
------------------------------------------------------------------------

Version Updates:

+ IgBLAST 1.9.0

Pipeline Changes:

+ Changed the default threshold detection method in ``shazam-threshold``
  to the smoothed density estimate with subsampling to 15000 sequences.
+ Fixed a bug wherein ``changeo-igblast`` was not reading the ``-b``
  argument.

Image Changes:

+ Added RDI R package.
+ Added CD-HIT.
+ Added AIRR python and R reference libaries.
+ Added git, BLAS, and LAPACK to base image.


Version 1.9.0:  April 22, 2018
------------------------------------------------------------------------

Version Updates:

+ alakazam 0.2.10
+ shazam 0.1.9

Pipeline Changes:

+ Added ``-l <model>`` argument to ``shazam-threshold`` to allow
  specification of the mixture model distributions to
  ``shazam::findThreshold``.

Image Changes:

+ Set Rcpp version for R package builds to ``0.12.16`` (from ``0.12.12``).


Version 1.8.0:  March 22, 2018
------------------------------------------------------------------------

Version Updates:

+ alakazam 0.2.9
+ changeo 0.3.12
+ presto 0.5.7

Pipeline Changes:

+ Removed an intermediate file and the ParseHeaders-rename step in
  ``presto-abseq``.
+ Modifed ``tigger-genotype`` to work with upcoming release of
  tigger v0.2.12.
+ Fixed parsing of output directory argument (``-o``) in
  ``preprocess-phix`` and ``changeo-clone``.

Image Changes:

+ Added sudo access for the magus (default) user.


Version 1.7.0:  February 6, 2018
------------------------------------------------------------------------

Version Updates:

+ changeo 0.3.11


Version 1.6.0:  January 29, 2018
------------------------------------------------------------------------

Version Updates:

+ prestor 0.0.4


Version 1.5.0:  January 17, 2018
------------------------------------------------------------------------

Version Updates:

+ presto 0.5.6


Version 1.4.0:  December 29, 2017
------------------------------------------------------------------------

Version Updates:

+ presto 0.5.5
+ phylip 3.697

Pipeline Changes:

+ Fixed a bug in ``presto-abseq`` preventing relative file paths from
  working with the ``-r`` argument.
+ ``changeo-igblast`` no longer terminates upon IgBLAST warnings.

Accessory Script Changes:

+ Fixed an output directory bug in ``fastq2fasta.py``.

Image Changes:

+ Added Stern, Yaari and Vander Heiden, et al 2014 primer sets.


Version 1.3.0:  October 17, 2017
------------------------------------------------------------------------

Version Updates:

+ changeo 0.3.9

Pipeline Changes:

+ Fixed a bug in ``presto-abseq`` preventing relative file paths from
  working with the ``-r`` argument.


Version 1.2.0:  October 05, 2017
------------------------------------------------------------------------

Version Updates:

+ changeo 0.3.8


Version 1.1.0:  September 22, 2017
------------------------------------------------------------------------

Version Updates:

+ alakazam 0.2.8
+ tigger 0.2.11
+ prestor 0.0.3

Image Changes:

+ Added ``preprocess-phix`` script that removes PhiX reads.
+ Added ``fetch_phix.sh`` script that downloads the PhiX174 genome.
+ Added ``builds`` script to record and report image build date and
  package changesets.
+ Added ``-x <coordinate system>`` argument to presto-abseq.
+ Forced install of Rcpp to be fixed at version 0.12.12.
+ Added ``/oasis`` mount point


Version 1.0.0:  August 08, 2017
------------------------------------------------------------------------

+ Initial meta-versioned image.