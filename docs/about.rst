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
    .. image:: _static/overview.svg


..
    .. graphviz::

        graph {
            graph[rankdir="LR"];

            a -- b;
            b -- c;
            c -- d;

            d -- e;
            d -- f;
            d -- f;
            d -- e;
            d -- e;
            d -- e;

            b -- d;
            b -- e;

            a [label="Data"];
            b [label="Preprocessing"];
            c [label="V(D)J Alignment"];
            d [label="Standardization"];
            e [label="Clonal Assignment"];
            f [label="Lineage Reconstruction"];
            g [label="Diversity"];
            h [label="Mutation Profiling"];
            i [label="Selection"];
            j [label="Genotyping"];

            k [label="pRESTO", href="http://presto.readthedocs.org"];
            l [label="Change-O", href="http://changeo.readthedocs.org"];
            m [label="Alakazam", href="http://kleinstein.bitbucket.org/alakazam"];
            n [label="Shazam", href="http://kleinstein.bitbucket.org/shazam"];
            o [label="TIgGER", href="http://clip.med.yale.edu/tigger"];
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