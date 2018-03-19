Docker Container Release Notes
========================================================================

Version 1.8.0:  March 19, 2018
------------------------------------------------------------------------

Version Updates:

+ changeo 0.3.12
+ presto 0.5.7

Pipeline Changes:

+ Removed an intermediate file and the ParseHeaders-rename step in
  ``presto-abseq``.
+ Added support for either ``.tsv`` or ``.tab`` extension to
  ``changeo-igblast`` and ``changeo-clone``.
+ Modifed ``tigger-genotype`` to work with TIgGER v0.2.12.
+ Fixed parsing of output directory argument (``-o``) in
  ``preprocess-phix`` and ``changeo-clone``.

Image Changes:

+ Added sudo access for the magus (default) user.


Version 1.7.0:  February 6, 2018
------------------------------------------------------------------------

Version Updates:

+ changeo 0.3.11


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

Image Changes:

+ Added Stern, Yaari and Vander Heiden, et al 2014 primer sets.
+ Fixed an output directory bug in ``fastq2fasta.py``.


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