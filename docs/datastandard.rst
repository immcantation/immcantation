Data Standards
===========================================================================================

Immcantation supports both the original Change-O standard and the new
`Adaptive Immune Receptor Repertoire (AIRR) <https://docs.airr-community.org/en/latest/index.html>`__
standard developed by the AIRR Community (`AIRR-C <https://www.antibodysociety.org/the-airr-community/>`__).
Both standards use tab-delimited file format files with a set of predefined column names.

AIRR Community Data Standard
-------------------------------------------------------------------------------------------

As of release 4.0.0, the default file format is the AIRR-C format, as described by the Rearrangement 
Schema (version 1.2). To learn about this format, the valid field names and their expected values, visit the 
AIRR-C `Rearrangement Schema documentation site <https://docs.airr-community.org/en/v1.2.1/datarep/rearrangements.html>`__.

    **Vander Heiden et al**
    AIRR Community Standardized Representations for Annotated Immune Repertoires.
    *Frontiers in Immunology.    9, 2206 (2018).*
    `doi\:10.3389/fimmu.2018.02206 <https://doi.org/10.3389/fimmu.2018.02206>`__
    
Change-O Format
-------------------------------------------------------------------------------------------

The Change-O format is the original common data format developed to enable the integration of 
the multiple tools ind the Immcantation framework. It is described in detail in the documentation
of `Change-O <https://changeo.readthedocs.io/en/latest/standard.html>`__.

    **Gupta NT\*, Vander Heiden JA\*, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH.**
    Change-O\: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data.
    *Bioinformatics 2015; doi\: 10.1093/bioinformatics/btv359*

Change of the default format in v4.0.0
--------------------------------------------------------------------------------------------

Release 4.0.0 introduces two main changes that can potentially break Immcantation pipelines. 
In this  section we explain these changes and provide an example to show how to update pipelines to
work with release 4.0.0.

The first chage is the adoption of the AIRR Standard as the default format expected by the 
tools, which implies default values in all functions and pipelines have been adjusted to use
AIRR Standard values. Users upgrading to 4.0.0, can experience that workflows that relied 
on default values now fail. The solution is to review the workflow and specify the correct 
values for the data format used.

.. parsed-literal::

   library(alakazam)
   library(shazam)
    
   # alakazam provides an example dataset in Change-O format
   db <- ExampleDbChangeo
   
   # Inspect the column names 
   > colnames(db)
   [1] "SEQUENCE_ID"          "SEQUENCE_IMGT"        "GERMLINE_IMGT_D_MASK"
   [4] "V_CALL"               "V_CALL_GENOTYPED"     "D_CALL"              
   [7] "J_CALL"               "JUNCTION"             "JUNCTION_LENGTH"     
   [10] "NP1_LENGTH"           "NP2_LENGTH"           "SAMPLE"              
   [13] "ISOTYPE"              "DUPCOUNT"             "CLONE"  
   
   # As of release 4.0.0, the command below doesn't work if the input data is in 
   # Change-O format, because the default values in `distToNearest` are 
   # AIRR Standard values:
   #    sequence_column="junction"
   #    vCallColumn="v_call"
   #    jCallColumn="j_call"
   # These values don't match column names in `db`, therefore the command doesn't work
   db <- distToNearest(db)      
   
   # The solution is to specify the correct values: 
   db <- distToNearest(db, sequence_column="JUNCTION", vCallColumn="V_CALL", jCallColumn="J_CALL")      

   
   
   


