.. _workflow-configuration:

######################
Workflow configuration
######################

HRIBO allows different customization to be able to handle different types of input data and customize the analysis.
On this page we explain the different options that can be set to easily customize the workflow.


Default workflow
================

We provide a default ``config.yaml`` with generally well-defined default values in the template folder of the HRIBO GitHub repository.


Biology Settings
================

This section contains general settings for the workflow.

adapter sequence
****************

.. code-block:: bash

    adapter: ""

The adapter sequence that will be removed using ``cutadapt``. It is important to add it if you know that your data was not trimmed before.

genome file
***********

.. code-block:: bash

    genome: "genome.fa"

The path to the genome file in fasta format. Only fasta format is supported. Make sure that the genome file has the same identifiers as the reference annotation file.

annotation file
***************

.. code-block:: bash

    annotation: "annotation.gff"

The path to the annotation file in gff format. Only gff format is supported. Make sure that the annotation file has the same identifiers as the reference genome file.

sample sheet
************

.. code-block:: bash

    samples: "HRIBO/samples.tsv"

The path to the sample sheet (example provided in the templates folder). The sample sheet describes the relation between the different samples used for the experiment.

alternative start codons
************************

.. code-block:: bash

    alternativestartcodons: ["GTG","TTG"]

A list of alternative start codons that will be used to predict ORFs. The default is ``["GTG","TTG"]``.


Differential Expression Settings
================================

differentialexpression
**********************

.. code-block:: bash

    differentialexpression: "on"

This option allows you to turn on or off differential expression analysis. If you do not have multiple conditions defined in the sample sheet and differential expression is activated, you will receive an error message.
Options are "on / off".

features
********

.. code-block:: bash

    features: ["CDS", "sRNA"]

This option allows you to specify which features should be used for differential expression analysis. Any feature that appears in your reference annotation is allowed.
We suggest using CDS and sRNA features.

contrasts
*********

.. code-block:: bash

    contrasts: ["treated1-untreated1", "treated2-untreated2"]

This option allows you to specify which contrasts should be used for differential expression analysis.
The order will affect the directionality of the log2FC values in the output files.

adjusted pvalue cutoff
**********************

.. code-block:: bash

    padjCutoff: 0.05

This option allows you to specify the adjusted pvalue cutoff for differential expression analysis. The default is 0.05.
All results will be present in the output, this will is soley used to add additional pre-filtered list in the output excel tables.

log2 fold change cutoff
***********************

.. code-block:: bash

    log2fcCutoff: 1.0

This option allows you to specify the log2 fold change cutoff for differential expression analysis. The default is 1.0.
All results will be present in the output, this will is soley used to add additional pre-filtered list in the output excel tables.
Only positive values are allowed. Respective negative values to determine down-regulation will be generated automatically by multiplying by -1.

ORF predictions
===============

.. code-block:: bash

    deepribo: "on"

Activating DeepRibo predictions will give you a different file with ORF predictions.
By experience, the top DeepRibo results tend to be better than those of reparation.
For archea, where reparation performs very poorly, DeepRibo is the preferred option.

.. warning:: DeepRibo cannot cope with genomes containing special ``IUPAC symbols``, ensure that your genome file contains only ``A``, ``G``, ``C``, ``T``, ``N`` symbols.


Read statistics Settings
========================

.. code-block:: bash

    readLengths: "10-80"

This option allows you to specify the read lengths that should be used for read statistics analysis. The default is ``10-80``.
We allow combinations of intervals and single values, e.g. ``10-40,55,70``.


Metagene Profiling Settings
===========================

There exist multiple options to customize the metagene profiling analysis.

Positions outside of the ORF
****************************

.. code-block:: bash

    positionsOutsideORF: 50

This option allows you to specify the number of positions outside of the ORF that should be considered for metagene profiling analysis. The default is ``50``.

Positions inside of the ORF
**************************

.. code-block:: bash

    positionsInsideORF: 150

This option allows you to specify the number of positions inside of the ORF that should be considered for metagene profiling analysis. The default is ``150``.
Genes that are shorter than the specified number of positions inside of the ORF will be ignored.

Filtering Methods
*****************

.. code-block:: bash

    filteringMethods: ["overlap", "length", "rpkm"]

This option allows you to specify which filtering methods should be used for metagene profiling analysis. The default is ``["overlap", "length", "rpkm"]``.

These methods are used to filter out genes that would cause artifacts in the metagene profiling analysis.

* The ``overlap`` method filters out genes that overlap with other genes.
* The ``length`` method filters out genes that are shorter than the threshold.
* The ``rpkm`` method filters out genes below the rpkm threshold.

Neighboring Genes
*****************

.. code-block:: bash

    neighboringGenes: 1

Paired-end support
==================

Until full paired-end support is added, we allow paired-end data by merging it into single-end data.
Therefore, we use the tool ``flash2`` to convert paired-end data to single-end data.

In order to use paired-end data, simply replace the ``Snakefile`` with the ``Snakefile_pairedend``.
This will now require a special ``samples_pairedend.tsv``, which is also available in the HRIBO templates folder.

+-----------+-----------+-----------+----------------------------+----------------------------+
|   method  | condition | replicate | fastqFile                  | fastqFile2                 |
+===========+===========+===========+============================+============================+
| RIBO      |  A        | 1         | fastq/RIBO-A-1_R1.fastq.gz | fastq/RIBO-A-1_R2.fastq.gz |
+-----------+-----------+-----------+----------------------------+----------------------------+
| RIBO      |  A        | 2         | fastq/RIBO-A-2_R1.fastq.gz | fastq/RIBO-A-2_R2.fastq.gz |
+-----------+-----------+-----------+----------------------------+----------------------------+
| RIBO      |  B        | 1         | fastq/RIBO-B-1_R1.fastq.gz | fastq/RIBO-B-1_R2.fastq.gz |
+-----------+-----------+-----------+----------------------------+----------------------------+
| RIBO      |  B        | 2         | fastq/RIBO-B-2_R1.fastq.gz | fastq/RIBO-B-2_R2.fastq.gz |
+-----------+-----------+-----------+----------------------------+----------------------------+
| RNA       |  A        | 1         | fastq/RNA-A-1_R1.fastq.gz  | fastq/RNA-A-1_R2.fastq.gz  |
+-----------+-----------+-----------+----------------------------+----------------------------+
| RNA       |  A        | 2         | fastq/RNA-A-2_R1.fastq.gz  | fastq/RNA-A-2_R2.fastq.gz  |
+-----------+-----------+-----------+----------------------------+----------------------------+
| RNA       |  B        | 1         | fastq/RNA-B-1_R1.fastq.gz  | fastq/RNA-A-1_R2.fastq.gz  |
+-----------+-----------+-----------+----------------------------+----------------------------+
| RNA       |  B        | 2         | fastq/RNA-B-2_R1.fastq.gz  | fastq/RNA-A-1_R2.fastq.gz  |
+-----------+-----------+-----------+----------------------------+----------------------------+

