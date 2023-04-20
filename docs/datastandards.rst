Data Standards
===========================================================================================

Immcantation supports both the original `Change-O standard <https://changeo.readthedocs.io/en/latest/standard.html>`__ and the new
`Adaptive Immune Receptor Repertoire (AIRR) <https://docs.airr-community.org/en/stable/>`__
standard developed by the AIRR Community (`AIRR-C <https://www.antibodysociety.org/the-airr-community/>`__).
Both standards use tab-delimited file formats with sets of specific predefined column names.

Change-O Standard
-------------------------------------------------------------------------------------------

The Change-O format is the original data format developed to enable the integration of
multiple tools in the Immcantation framework. It is described in detail (along with the corresponding AIRR-C Standard equivalents) in the
`Change-O package documentation <https://changeo.readthedocs.io/en/latest/standard.html>`__.

    **Gupta NT\*, Vander Heiden JA\*, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH.**
    Change-O\: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data.
    *Bioinformatics 2015*; `doi\: 10.1093/bioinformatics/btv359 <https://doi.org/10.1093/bioinformatics/btv359>`__

AIRR Community Standard
-------------------------------------------------------------------------------------------

The default file format for all functions in Immcantation is the AIRR-C format as of release 4.0.0.
To learn more about this format (including the valid field names and their expected values), visit the
`AIRR-C Rearrangement Schema documentation <https://docs.airr-community.org/en/stable/datarep/rearrangements.html>`__.
The `Change-O package documentation <https://changeo.readthedocs.io/en/latest/standard.html>`__ contains a table with mappings
between both standards. Some of the most frequently used translations are:

+------------------------+----------------+
| AIRR                   | Change-O       |
+========================+================+
| sequence_id            | SEQUENCE_ID    |
+------------------------+----------------+
| sequence               | SEQUENCE_INPUT |
+------------------------+----------------+
| sequence_alignment     | SEQUENCE_IMGT  |
+------------------------+----------------+
| productive             | FUNCTIONAL     |
+------------------------+----------------+
| v_call                 | V_CALL         |
+------------------------+----------------+
| d_call                 | D_CALL         |
+------------------------+----------------+
| j_call                 | J_CALL         |
+------------------------+----------------+
| junction_length        | JUNCTION_LENGTH|
+------------------------+----------------+
| junction               | JUNCTION       |
+------------------------+----------------+
| germline_alignment     | GERMLINE_IMGT  |
+------------------------+----------------+
| clone_id               | CLONE          |
+------------------------+----------------+

    **Vander Heiden et al**
    AIRR Community Standardized Representations for Annotated Immune Repertoires.
    *Frontiers in Immunology. 9, 2206 (2018).*
    `doi\:10.3389/fimmu.2018.02206 <https://doi.org/10.3389/fimmu.2018.02206>`__


Potential Workflow Changes
--------------------------------------------------------------------------------------------

Release 4.0.0 introduces two main changes that can potentially break existing Immcantation workflows.
In this section, we explain these changes, give solutions and provide an example to show how to update workflows in order to properly
work with release 4.0.0.

The first change is the adoption of the AIRR Standard as the default format expected by the
tools (note that Change-O is still available as an option). The default values in **all functions and pipelines**
have been adjusted to use this standard. Users upgrading to 4.0.0 may find that workflows that relied upon
default values now fail. *The solution is to review the workflow and specify the correct values for the data format being used.*

The second change that can break workflows is that **all outputs now use lowercase column names** for style consistency with the AIRR Standard
format. This means that user workflows that expect columns to be in uppercase will now
break. *The solution is to update the code to use the current lowercase values.*

The following R-based example demonstrates how to fix broken workflows as a result of these two changes:

.. code-block:: R

    > library(alakazam)
    > library(shazam)

    # alakazam provides an example dataset in Change-O format
    > db <- ExampleDbChangeo

    # Inspect the column names
    > colnames(db)
    [1] "SEQUENCE_ID"          "SEQUENCE_IMGT"        "GERMLINE_IMGT_D_MASK"
    [4] "V_CALL"               "V_CALL_GENOTYPED"     "D_CALL"
    [7] "J_CALL"               "JUNCTION"             "JUNCTION_LENGTH"
    [10] "NP1_LENGTH"           "NP2_LENGTH"           "SAMPLE"
    [13] "ISOTYPE"              "DUPCOUNT"             "CLONE"

    # CHANGE 1: default values follow the AIRR Standard specification
    > db <- distToNearest(db)
    Error in distToNearest(db) : The column junction was not found

    # As of release 4.0.0, the `distToNearest` command above doesn't work if the input data
    # is in Change-O format because the default values are now AIRR Standard values:
    #    sequenceColumn="junction"
    #    vCallColumn="v_call"
    #    jCallColumn="j_call"
    # These values don't match the column names in `db` as previously seen, so the command doesn't work

    # The solution is to specify the actual column names:
    > db <- distToNearest(db, sequenceColumn="JUNCTION",
                            vCallColumn="V_CALL",
                            jCallColumn="J_CALL")
    > colnames(db)
    [1] "SEQUENCE_ID"          "SEQUENCE_IMGT"        "GERMLINE_IMGT_D_MASK"
    [4] "V_CALL"               "V_CALL_GENOTYPED"     "D_CALL"
    [7] "J_CALL"               "JUNCTION"             "JUNCTION_LENGTH"
    [10] "NP1_LENGTH"           "NP2_LENGTH"           "SAMPLE"
    [13] "ISOTYPE"              "DUPCOUNT"             "CLONE"
    [16] "dist_nearest"

    # CHANGE 2: outputs are generated using lower case
    > threshold <- findThreshold(db$DIST_NEAREST)
    Error in h.ucv.default(unique(distances), 4) :
      argument 'x' must be numeric and need at least 3 data points
    In addition: Warning message:
    Unknown or uninitialized column: 'DIST_NEAREST'.

    # In previous releases, `distToNearest` added the column `DIST_NEAREST` to `db`.
    # As of release 4.0.0, it adds `dist_nearest`, so the command above
    # doesn't work, because `db` doesn't have a column named `DIST_NEAREST`

    # The solution is to update the function call to use the correct name:
    > threshold <- findThreshold(db$dist_nearest)


Convert between Change-O and AIRR-C format
-------------------------------------------------------------------------------------------
The default file format is the AIRR-C format. However, Immcantation provides a script ConvertDb in Change-O package to convert the file from AIRR-C format to the legacy Change-O standard. For example, to convert a file named sample1_airr.tsv in AIRR-C format to Change-O format, you can run:

.. code-block:: R

    > ConvertDb changeo -d sample1_airr.tsv -o sample1_changeo.tab

The output file sample1_changeo.tab is in Change-O format.

In a similar way, you can also use ConvertDb to convert a Change-O file to an AIRR-C file:

.. code-block:: R

    > ConvertDb airr -d sample1_changeo.tab -o sample1_airr.tsv
