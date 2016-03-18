.. Ig Tools documentation master file
.. keywords:  rep-seq, immuno-seq, repertoire sequencing, BCR, TCR, Ig, IgOme, AIRR, adaptive immunity,
    somatic hypermutation, vdj-seq, AbSeq, AbPair, iRepertoire, VDJ, Immunoglobulin

.. toctree::
   :maxdepth: 2

   `pRESTO <presto.readthedocs.org>`_
   `Change-O <changeo.readthedocs.org>`_
   `aLAkazam <bitbuckert.org/kleinstein/alakazam>`_
   `SHazaM <bitbuckert.org/kleinstein/shazam>`_
   `TIgGER <bitbuckert.org/kleinstein/tigger>`_


Immcantation: An Integrated Framework for Adaptive Immune Receptor Repertoire Analysis
================================================================================

Advances in high-throughput sequencing technologies now allow for large-scale characterization of
B cell receptor (BCR) and T cell receptor (TCR) repertoires. The high germline and somatic diversity of the
adaptive immune receptor repertoire (AIRR) presents challenges for biologically meaningful analysis -
requiring the development of specialized computational methods.

Immcantation is a start-to-finish analytical ecosystem for high-throughput AIRR-seq datasets.
Beginning from raw reads, Python and R packages are provided for pre-processing, population structure determination,
and repertoire analysis.


Flowchart
--------------------------------------------------------------------------------

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

Links
--------------------------------------------------------------------------------

pRESTO: Pre-processing and UID/UMI analysis
Change-O: IMGT/HighV-QUEST and IgBLAST import, clonal grouping
TIgGER: Novel allele detection and genotyping
aLAkazam: Lineage analysis, diversity, amino acid properties
SHazaM: Somatic hypermutation analysis

How to cite
--------------------------------------------------------------------------------

pRESTO
Change-O

Additional citations for specific methods within aLAkazam, SHazaM and TIgGER may be determined using the ``citation()`` function.


.. Publications that use Immcantation:


Contact Info
--------------------------------------------------------------------------------

steve
website


