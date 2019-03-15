Release Notes
========================================================================

Version 2.8.0:  March 19, 2018
------------------------------------------------------------------------

Version Updates:

+ igblast 1.13.0


Version 2.7.0:  February 1, 2019
------------------------------------------------------------------------

Version Updates:

+ presto 0.5.11
+ changeo 0.4.5
+ shazam 0.1.11
+ blast 2.8.1


Version 2.6.0:  December 9, 2018
------------------------------------------------------------------------

Version Updates:

+ igblast 1.12.0

Pipeline Changes:

+ Added ``-i`` argument to ``changeo-igblast`` to allow retention of
  partial alignments.
  
Image Changes:

+ Base system changed to Fedora 29.
+ Moved setup of R package build environment to base image.


Version 2.5.0:  November 1, 2018
------------------------------------------------------------------------

Version Updates:

+ igblast 1.11.0
+ muscle 3.8.425
+ vsearch 2.9.1

Image Changes:

+ Added error checking to ``versions report`` command.


Version 2.4.0:  October 27, 2018
------------------------------------------------------------------------

Version Updates:

+ changeo 0.4.4


Version 2.3.0:  October 21, 2018
------------------------------------------------------------------------

Version Updates:

+ presto 0.5.10
+ changeo 0.4.3
+ tigger 0.3.1

Image Changes:

+ Added scoper R package.
+ Added IgPhyML.
+ Removed strict Rcpp version requirement (was fixed at ``0.12.16``).
+ Added libGL and libGLU to base image.


Version 2.2.0:  October 5, 2018
------------------------------------------------------------------------

Version Updates:

+ tigger 0.3.0
+ airr python library 1.2.1

Pipeline Changes:

+ Fixed compression error messages in ``changeo-igblast`` and
  ``changeo-clone``.
+ Removed support for tigger versions below 0.3.0 from
  ``tigger-genotype``.

Image Changes:

+ Adjusted version/changeset detection and output in the
  ``versions report`` and ``builds report`` commands.


Version 2.1.0:  September 20, 2018
------------------------------------------------------------------------

Version Updates:

+ alakazam 0.2.11
+ shazam 0.1.10
+ prestor 0.0.5
+ vsearch 2.8.4
+ BLAST 2.7.1
+ IgBLAST 1.10.0

Pipeline Changes:

+ Subsampling is no longer performed by default in ``shazam-threshold``.

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
