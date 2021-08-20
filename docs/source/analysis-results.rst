.. _analysis-results:

#####################
Analysis result files
#####################

The important files in this workflow are listed and explained below.

ORF Predictions
===============

The output files containing information about predicted Open Reading Frames, these also contain novel predictions.

predictions_reparation.xlsx
***************************

This file contains all ``reparation`` ORF predictions.

+-------------------------------------------+-----------------------------------------------------------------------------+
| Column name                               | Description                                                                 |
+===========================================+=============================================================================+
| Identifier                                | Unique identifier describing the entry.                                     |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Genome                                    | The genome accession identifier.                                            |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Source                                    | The source of the ORF. (here reparation)                                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Feature                                   | The feature of the ORF (here CDS)                                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start                                     | The start position of the ORF.                                              |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop                                      | The stop position of the ORF.                                               |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Strand                                    | The strand of the ORF. (+/-)                                                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Pred_probability                          | The probability of this ORF (0.5-1, 1 being the best value)                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Locus_tag                                 | If the detected ORF is already in the annotation, this gives its locus tag. |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Old_locus_tag                             | The old locus tag of a gene (if available in the annotation)                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Name                                      | If the detected ORF is already in the annotation, this gives its name.      |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Length                                    | The length of the ORF.                                                      |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Codon_count                               | The number of codons in the ORF. (length/3)                                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_TE       | The translational efficiency for the given sample.                          |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_rpkm     | The RPKM for the given sample.                                              |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Evidence                                  | The <condition>-<replicate> sample in which this ORF was predicted.         |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start_codon                               | The start codon of the ORF.                                                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop_codon                                | The stop codon of the ORF.                                                  |
+-------------------------------------------+-----------------------------------------------------------------------------+
| 15nt_upstream                             | The 15nt upstream of the start codon                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Nucleotide_seq                            | The nucleotide sequence of the ORF.                                         |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Aminoacid_seq                             | The amino acid sequence of the ORF.                                         |
+-------------------------------------------+-----------------------------------------------------------------------------+

predictions_reparation.gff
**************************

An annotation file in ``.gff3`` format containing all predictions of ``reparation`` for visualization in a genome browser.


predictions_deepribo.xlsx
*************************

.. note:: These files are only available when activating DeepRibo predictions in the ``config.yaml``. (see :ref:`workflow-configuration <source/workflow-configuration:workflow-configuration`>)

This file contains all ``DeepRibo`` ORF predictions.

+-------------------------------------------+---------------------------------------------------------------------------------+
| Column name                               | Description                                                                     |
+===========================================+=================================================================================+
| Identifier                                | Unique identifier describing the entry.                                         |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Genome                                    | The genome accession identifier.                                                |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Source                                    | The source of the ORF. (here reparation)                                        |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Feature                                   | The feature of the ORF (here CDS)                                               |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Start                                     | The start position of the ORF.                                                  |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Stop                                      | The stop position of the ORF.                                                   |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Strand                                    | The strand of the ORF. (+/-)                                                    |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Pred_value                                | The value DeepRibo attributes the given prediction.                             |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Pred_rank                                 | The rank calculated from the prediction value. (the best prediction has rank 1) |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Novel_rank                                | A special ranking involving only novel ORFs that are not in the annotation.     |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Locus_tag                                 | If the detected ORF is already in the annotation, this gives its locus tag.     |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Old_locus_tag                             | The old locus tag of a gene (if available in the annotation)                    |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Name                                      | If the detected ORF is already in the annotation, this gives its name.          |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Length                                    | The length of the ORF.                                                          |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Codon_count                               | The number of codons in the ORF. (length/3)                                     |
+-------------------------------------------+---------------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_TE       | The translational efficiency for the given sample.                              |
+-------------------------------------------+---------------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_rpkm     | The RPKM for the given sample.                                                  |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Evidence                                  | The <condition>-<replicate> sample in which this ORF was predicted.             |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Start_codon                               | The start codon of the ORF.                                                     |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Stop_codon                                | The stop codon of the ORF.                                                      |
+-------------------------------------------+---------------------------------------------------------------------------------+
| 15nt_upstream                             | The 15nt upstream of the start codon                                            |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Nucleotide_seq                            | The nucleotide sequence of the ORF.                                             |
+-------------------------------------------+---------------------------------------------------------------------------------+
| Aminoacid_seq                             | The amino acid sequence of the ORF.                                             |
+-------------------------------------------+---------------------------------------------------------------------------------+

predictions_deepribo.gff
************************

.. note:: These files are only available when activating DeepRibo predictions in the ``config.yaml``. (see :ref:`workflow-configuration <source/workflow-configuration:Workflow configuration`>)

An annotation file in ``.gff3`` format containing all predictions of *DeepRibo* for visualization in a genome browser.


Quality control
===============

This comprises all files that can help to perform quality control on all input samples.

multiqc_report.html
*******************

The multiQC report collects information from different tools, including ``fastQC`` and ``subread featurecounts``.
The general statistics give an overview over:

*	the number of duplicates
*	the GC content
*	the average read lengths
*	the number of reads (in millions)

These statistics are collected after each processing step of our pipeline.

*	**raw:** the unprocessed data
*	**trimmed:** the data after trimming the adapter sequences
*	**mapped:** the data after mapping with Segemehl
*	**unique:** the data after removing multi-mapping reads
*	**norRNA:** the data after filtering out the rRNA

Further, feature counts are provided for different features from the annotation file. (i.e. how many reads map to each feature)
This includes, all(featurecount), rRNA, norRNA(after filtering), tRNA and ncRNA.
Following is a fastQC report including sequence counts, sequence quality histograms, per sequence quality scores, per base sequence content, per sequence GC content, per base N content, sequence length distribution, sequence duplication levels, overrepresented features, adapter content and a status overview.


heatmap_SpearmanCorr_readCounts.pdf
***********************************

Spearman correlation coefficients of read counts. The dendrogram indicates which samples read counts are most similar to each other.
Since there should be always a higher correlation between experiments with the same condition and experiment type (e.g. replicates) and not others, this is a rapid way to quality-control the labeling/consistency of input data.

annotation_total.xlsx
*********************

This file contains detailed measures for every feature in the input annotation using read counts including multi-mapping reads.

+-------------------------------------------+-----------------------------------------------------------------------------+
| Column name                               | Description                                                                 |
+===========================================+=============================================================================+
| Identifier                                | Unique identifier describing the entry.                                     |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Genome                                    | The genome accession identifier.                                            |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Source                                    | The source of the annotated feature.                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Feature                                   | The feature of the annotated feature.                                       |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start                                     | The start position of the annotated feature.                                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop                                      | The stop position of the annotated feature.                                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Strand                                    | The strand of the annotated feature. (+/-)                                  |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Locus_tag                                 | The locus tag of the annotated feature. (if available)                      |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Old_locus_tag                             | The old locus tag of a gene (if available in the annotation)                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Name                                      | The name of the annotated feature. (if available)                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Length                                    | The length of the annotated feature.                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Codon_count                               | The number of codons in the annotated feature. (length / 3)                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_TE       | The translational efficiency for the given sample.                          |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_rpkm     | The RPKM for the given sample. (ReadsPerKilobaseMillion)                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start_codon                               | The start codon of the annotated feature.                                   |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop_codon                                | The stop codon of the annotated feature.                                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| 15nt_upstream                             | The 15nt upstream of the start codon                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Nucleotide_seq                            | The nucleotide sequence of the annotated feature.                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Aminoacid_seq                             | The amino acid sequence of the annotated feature.                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Product                                   | The product of the annotated feature. (if available)                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Note                                      | The note of the annotated feature. (if available)                           |
+-------------------------------------------+-----------------------------------------------------------------------------+

total_read_counts.xlsx
**********************

This file shows the overall read-counts for each feature annotated in the user-provided annotation, after mapping and before removal of multi-mapping reads.

annotation_unique.xlsx
**********************

This file contains detailed measures for every feature in the input annotation using read counts after removal of multi-mapping reads.

+-------------------------------------------+-----------------------------------------------------------------------------+
| Column name                               | Description                                                                 |
+===========================================+=============================================================================+
| Identifier                                | Unique identifier describing the entry.                                     |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Genome                                    | The genome accession identifier.                                            |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Source                                    | The source of the annotated feature.                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Feature                                   | The feature of the annotated feature.                                       |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start                                     | The start position of the annotated feature.                                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop                                      | The stop position of the annotated feature.                                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Strand                                    | The strand of the annotated feature. (+/-)                                  |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Locus_tag                                 | The locus tag of the annotated feature. (if available)                      |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Old_locus_tag                             | The old locus tag of a gene (if available in the annotation)                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Name                                      | The name of the annotated feature. (if available)                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Length                                    | The length of the annotated feature.                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Codon_count                               | The number of codons in the annotated feature. (length / 3)                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_TE       | The translational efficiency for the given sample.                          |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_rpkm     | The RPKM for the given sample. (ReadsPerKilobaseMillion)                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start_codon                               | The start codon of the annotated feature.                                   |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop_codon                                | The stop codon of the annotated feature.                                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| 15nt_upstream                             | The 15nt upstream of the start codon                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Nucleotide_seq                            | The nucleotide sequence of the annotated feature.                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Aminoacid_seq                             | The amino acid sequence of the annotated feature.                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Product                                   | The product of the annotated feature. (if available)                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Note                                      | The note of the annotated feature. (if available)                           |
+-------------------------------------------+-----------------------------------------------------------------------------+

unique_read_counts.xlsx
***********************

This file shows the overall read-counts for each feature annotated in the user-provided annotation, after mapping and after removal of multi-mapping reads.


genome-browser
==============

The files that can be used for visualization in a genome browser.

updated_annotation.gff
**********************

A gff track containing both the original annotation together with the new predictions by reparation.

potentialStartCodons.gff
************************

A genome browser track with all possible start codons.

potentialStopCodons.gff
***********************

A genome browser track with all possible stop codons.

potentialRibosomeBindingSite.gff
********************************

A genome browser track with possible ribosome binding sites.

potentialAlternativeStartCodons.gff
***********************************

A genome browser track with alternative start codons.

BigWig coverage files
*********************

We offer many different single nucleotide mapping bigwig files for genome browser visualization.
These files are available for different regions and performed with different methods.

* **global:** full read is mapped
* **centered:** region around the center.
* **threeprime:** region around the three prime end.
* **fiveprime:** region around the five prime end.

These are all available with the following normalization methods:

* **raw:** raw, unprocessed files. This should only be used to check the coverage of a single file. It should not be used to compare to other files.
* **min:** normalized by number of minimal total reads per sample (factor = min. number of reads / number of reads). This is the recommended normalization when comparing different samples from the same experiment.
* **mil:** normalized by 1000000 (factor = 1000000 / number of reads). This is the recommended normalization when comparing different samples from the different experiments.

Differential Expression
=======================

Files related to the differential expression analysis.

riborex/<contrast>_sorted.xlsx
******************************

Table containing all differential expression results from *riborex*.

riborex/<contrast>_significant.xlsx
***********************************

Table containing significant differential expression results from *riborex* (pvalue < 0.05).

xtail/<contrast>_sorted.xlsx
****************************

Table containing all differential expression results from *xtail*.

xtail/<contrast>_significant.xlsx
*********************************

Table containing significant differential expression results from *xtail* (pvalue < 0.05).

xtail/r_<contrast>.pdf
**********************

This figure shows the RPF-to-mRNA ratios in two conditions, where the position
of each gene is determined by its RPF-to-mRNA ratio (log2R) in two conditions,
represented on the x-axis and y-axis respectively. The points will be color-coded with
the pvalue final obtained with xtail (more significant p values having darker color)

* **blue:** for genes with log2R larger in first condition than second condition.
* **red:** for genes with log2R larger in second condition than the first condition.
* **green:** for genes with log2R changing homodirectionally in two condition.
* **yellow:** for genes with log2R changing antidirectionally in two condition.

xtail/fc_<contrast>.pdf
***********************

This figure shows the result of the differential expression at the two expression levels,
where each gene is a dot whose position is determined by its log2 fold change (log2FC)
of transcriptional level (mRNA log2FC), represented on the x-axis, and the log2FC
of translational level (RPF log2FC), represented on the y-axis. The points will be
color-coded with the pvalue final obtained with xtail (more significant p values having
darker color)

* **blue:** for genes whos mRNA log2FC larger than 1 (transcriptional level).
* **red:** for genes whos RPF log2FC larger than 1 (translational level).
* **green:** for genes changing homodirectionally at both level.
* **yellow:** for genes changing antidirectionally at two levels.

Metagene Profiling Analysis
===========================

Please refer to the :ref:`metagene profiling <source/metagene-profiling:Metagene profiling>` page for further details.

<accession>_Z.Y_profiling.xlsx/tsv
**********************************

The table shows for a range of specific read lengths, how many reads on average over all start codons
in the genome have been mapped per nucleotide. The nucleotides range from 100 nucleotides upstream
of the start codon to 399 nucleotides downstream. The read counts are either raw or normalized by average read count per nucleotide, for the range around the start codon. Moreover different single nucleotide mapping variants are considered,
where only the 5', 3' or centered region of the read is counted.


<accession>_Z.Y_profiling.pdf
*****************************


Additional output
=================

samples.xlsx
************

An excel representation of the input sample file.

manual.pdf
**********

A PDF format file giving some explanations about the output files, contained in the final result report.

overview.xlsx
*************

An overview table containing all information gathered from the prediction tools and differential expression analysis.
The contents of this table change depending on which :ref:`options <source/workflow-configuration:Workflow configuration>` are set.
The overview table for the default workflow will contain annotation. reparation, deepribo and differential expression output.

+-------------------------------------------+-----------------------------------------------------------------------------+
| Column name                               | Description                                                                 |
+===========================================+=============================================================================+
| Identifier                                | Unique identifier describing the entry.                                     |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Genome                                    | The genome accession identifier.                                            |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start                                     | The start position of the ORF.                                              |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop                                      | The stop position of the ORF.                                               |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Strand                                    | The strand of the ORF. (+/-)                                                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Locus_tag                                 | The locus tag of ORF. (if not novel)                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Overlapping_genes                         | Genes that overlap with the predicted ORF                                   |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Old_locus_tag                             | The old locus tag of a gene (if available in the annotation)                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Name                                      | The name of the ORF. (if not novel)                                         |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Gene_name                                 | The name of the ORFs associated gene feature. (if not novel)                |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Length                                    | The length of the ORF.                                                      |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Codon_count                               | The number of codons in the ORF. (length / 3)                               |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Start_codon                               | The start codon of the annotated feature.                                   |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Stop_codon                                | The stop codon of the annotated feature.                                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| 15nt_upstream                             | The 15nt upstream of the start codon                                        |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Nucleotide_seq                            | The nucleotide sequence of the annotated feature.                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Aminoacid_seq                             | The amino acid sequence of the annotated feature.                           |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_TE       | The translational efficiency for the given sample.                          |
+-------------------------------------------+-----------------------------------------------------------------------------+
| <method>-<condition>-<replicate>_rpkm     | The RPKM for the given sample. (ReadsPerKilobaseMillion)                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Evidence_reparation                       | The sample this ORF was predicted in (for reparation)                       |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Reparation_probability                    | The probability of this ORF (0.5-1, 1 being the best value)                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Evidence_deepribo                         | The sample this ORF was predicted in (for deepribo)                         |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Deepribo_rank                             | The deepribo rank for this ORF. (1 being the best value, 999999 undefined)  |
+-------------------------------------------+-----------------------------------------------------------------------------+
| Deepribo_score                            | The score the deepribo rank is based on.                                    |
+-------------------------------------------+-----------------------------------------------------------------------------+
| riborex_pvalue                            | The pvalue (determined by riborex)                                          |
+-------------------------------------------+-----------------------------------------------------------------------------+
| riborex_pvalue_adjusted                   | The adjusted pvalue (determined by riborex)                                 |
+-------------------------------------------+-----------------------------------------------------------------------------+
| riborex_log2FC                            | The log2FC (determined by riborex)                                          |
+-------------------------------------------+-----------------------------------------------------------------------------+
| xtail_pvalue                              | The pvalue (determined by xtail)                                            |
+-------------------------------------------+-----------------------------------------------------------------------------+
| xtail_pvalue_adjusted                     | The adjusted pvalue (determined by xtail)                                   |
+-------------------------------------------+-----------------------------------------------------------------------------+
| xtail_log2FC                              | The log2FC (determined by xtail)                                            |
+-------------------------------------------+-----------------------------------------------------------------------------+
