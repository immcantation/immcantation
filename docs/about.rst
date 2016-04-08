.. _Intro:

Overview
==========================================================================================

Advances in high-throughput sequencing technologies now allow for large-scale characterization of
B cell receptor (BCR) and T cell receptor (TCR) repertoires. The high germline and somatic diversity of the
adaptive immune receptor repertoire (AIRR) presents challenges for biologically meaningful analysis -
requiring the development of specialized computational methods.

Immcantation is a start-to-finish analytical ecosystem for high-throughput AIRR-seq datasets.
Beginning from raw reads, Python and R packages are provided for pre-processing, population structure determination,
and repertoire analysis.


.. raw:: html

    <object data="_static/overview.svg" type="image/svg+xml"></object>

..
    .. graphviz::

        graph {
            a -- b;
            b -- c;
            b -- f;
            f -- c;
            c -- d;
            c -- e;

            a [shape="rectangle", label="Read Assembly & Quality Control", href="http://presto.readthedocs.org", target="_top", image="_static/green_button.png"];
            b [shape="rectangle", label="V(D)J Alignment"];
            c [shape="rectangle", label="Clonal Clustering", href="http://changeo.readthedocs.org", target="_top"];
            d [shape="rectangle", label="Alakazam", href="http://kleinstein.bitbucket.org/alakazam", target="_top"];
            e [shape="rectangle", label="Shazam", href="http://kleinstein.bitbucket.org/shazam", target="_top"];
            f [shape="rectangle", label="Genotyping", href="http://clip.med.yale.edu/tigger", target="_top"];
        }

.. _Cite:

How to Cite
==========================================================================================

To cite the pRESTO software package in publications please use:

    **Vander Heiden JA\*, Yaari G\*, Uduman M, Stern JNH, O'Connor KC, Hafler DA, Vigneault F, Kleinstein SH.**
    pRESTO\: a toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires.
    *Bioinformatics 2014; doi\: 10.1093/bioinformatics/btu138*

To cite the Change-O, Alakazam, SHazaM and TIgGER software package in publications please use:

    **Gupta NT\*, Vander Heiden JA\*, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH.**
    Change-O\: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data.
    *Bioinformatics 2015; doi\: 10.1093/bioinformatics/btv359*

Additional citations for specific methods within Alakazam, SHazaM and TIgGER may be determined
using the ``citation()`` function witin R.

.. _Contact:

Contact Information
==========================================================================================

If you have questions you can email
`Steven Kleinstein <mailto:steven.kleinstein@yale.edu>`__.

For additional computational immunology software from the Kleinstein Lab see our
`website <http://medicine.yale.edu/lab/kleinstein/software/>`__.

.. Publications that use Immcantation: