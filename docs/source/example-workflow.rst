.. _example-workflow:

################
Example workflow
################

The retrieval of input files and running the workflow locally and on a server cluster via a queuing system is demonstrated using an example with data available from `NCBI  <https://www.ncbi.nlm.nih.gov/>`_.
This dataset is available under the accession number *PRJNA379630*.

.. note:: This dataset, featuring *Pseudomonas aeruginosa PAO1*, was chosen for this tutorial as it contains both multiple conditions and replicates.
.. note:: Ensure that you have **miniconda3** and **singularity** installed and a conda environment set-up. Please refer to the :ref:`overview <overview:Tools>` for details on the installation.

Setup
=====

First of all, we start by creating the project directory and changing to it.

.. code-block:: bash

    $ mkdir project
    $ cd project

We then download the latest version of HRIBO into the newly created project folder and unpack it.

.. code-block:: bash

   $ wget https://github.com/RickGelhausen/HRIBO/archive/1.2.0.tar.gz
   $ tar -xzf 1.2.0.tar.gz; mv HRIBO-1.2.0 HRIBO; rm 1.2.0.tar.gz;

Retrieve and prepare input files
================================

Before starting the workflow, we have to acquire and prepare several input files. These files are the annotation file, the genome file, the fastq files, the configuration file and the sample sheet.

Annotation and genome files
***************************
First, we want to retrieve the annotation file and the genome file. In this case, we can find both on `NCBI  <https://www.ncbi.nlm.nih.gov/genome/187?genome_assembly_id=299953>`_ using the accession number *NC_002516.2*.

.. image:: images/todo.png
    :scale: 50%
    :align: center

On this page, we can directly retrieve both files by clicking on the according download links next to the file descriptions. Alternatively, you can directly download them using the following commands:

.. code block:: bash

    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.gff.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz

Then, we unpack and rename both files.

.. code block:: bash

    gunzip GCF_000006765.1_ASM676v1_genomic.gff.gz && mv GCF_000006765.1_ASM676v1_genomic.gff annotation.gff
    gunzip GCF_000006765.1_ASM676v1_genomic.fna.gz && mv GCF_000006765.1_ASM676v1_genomic.fna genome.fa

.. note:: We chose the default names for both annotation and genome file. You can also specify another path to the annotation or genome file in the config.yaml **TODO link**.

.fastq files
************

Next, we want to acquire the fastq files. The fastq files are available under the accession number *PRJNA379630* on `NCBI  <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA379630>`.
The files have to be downloaded using `The Sequence Read Archive (SRA)  <https://www.ncbi.nlm.nih.gov/sra/docs/>`.
There are multiple ways of downloading files from SRA as explained `here  <https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/>`.

As we already have conda installed, the easiest way is to install the *sra-tools*:

.. code-block:: bash

    conda create -n sra-tools -c bioconda -c conda-forge sra-tools

This will create a conda environment containing the sra-tools. Using these, we can simply pass the SRA identifiers and download the data:

.. code-block:: bash

    conda activate sra-tools;
    fasterq-dump SRR5356908; pigz -p 2 SRR5356908.fastq; mv SRR5356908.fastq.gz RNA-PAO1-gly-1.fastq.gz;
    fasterq-dump SRR5356906; pigz -p 2 SRR5356906.fastq; mv SRR5356906.fastq.gz RNA-PAO1-gly-2.fastq.gz;
    fasterq-dump SRR5356907; pigz -p 2 SRR5356907.fastq; mv SRR5356907.fastq.gz RIBO-PAO1-gly-1.fastq.gz;
    fasterq-dump SRR5356905; pigz -p 2 SRR5356905.fastq; mv SRR5356905.fastq.gz RIBO-PAO1-gly-2.fastq.gz;
    fasterq-dump SRR5356902; pigz -p 2 SRR5356902.fastq; mv SRR5356902.fastq.gz RNA-PAO1-nalk-1.fastq.gz;
    fasterq-dump SRR5356900; pigz -p 2 SRR5356900.fastq; mv SRR5356900.fastq.gz RNA-PAO1-nalk-2.fastq.gz;
    fasterq-dump SRR5356901; pigz -p 2 SRR5356901.fastq; mv SRR5356901.fastq.gz RIBO-PAO1-nalk-1.fastq.gz;
    fasterq-dump SRR5356899; pigz -p 2 SRR5356899.fastq; mv SRR5356899.fastq.gz RIBO-PAO1-nalk-2.fastq.gz;


.. note:: Due to the runtime of several tools, especially the mapping by segemehl. We only chose two conditions and two replicates from this dataset.

.. note:: If you have a bad internet connection, this step might take up to several hours. You might want to consider downloading only one replicate per condition, or use your own data instead.

This will download compressed files for each of the required *.fastq* files. We will move them into a folder called *fastq*.

.. code-block:: bash

    mv *.fastq.gz fastq


Sample sheet and configuration file
***********************************

Finally, we will prepare the configuration file (*config.yaml*) and the sample sheet (*samples.tsv*). We start by copying templates for both files from the *HRIBO/templates/* into the *HRIBO/* folder.

.. code-block:: bash

    $ cp HRIBO/templates/bam-samples.tsv HRIBO/
    $ mv HRIBO/bam-samples.tsv HRIBO/samples.tsv

The sample file looks as follows:

+-----------+-----------+-----------+-------------------------+
|   method  | condition | replicate | inputFile               |
+===========+===========+===========+=========================+
| RIBO      |  A        | 1         | fastq/RIBO-A-1.fastq.gz |
+-----------+-----------+-----------+-------------------------+
| RIBO      |  A        | 2         | fastq/RIBO-A-2.fastq.gz |
+-----------+-----------+-----------+-------------------------+
| RIBO      |  B        | 1         | fastq/RIBO-B-1.fastq.gz |
+-----------+-----------+-----------+-------------------------+
| RIBO      |  B        | 2         | fastq/RIBO-B-2.fastq.gz |
+-----------+-----------+-----------+-------------------------+
| RNA       |  A        | 1         | fastq/RNA-A-1.fastq.gz  |
+-----------+-----------+-----------+-------------------------+
| RNA       |  A        | 2         | fastq/RNA-A-2.fastq.gz  |
+-----------+-----------+-----------+-------------------------+
| RNA       |  B        | 1         | fastq/RNA-B-1.fastq.gz  |
+-----------+-----------+-----------+-------------------------+
| RNA       |  B        | 2         | fastq/RNA-B-2.fastq.gz  |
+-----------+-----------+-----------+-------------------------+

.. note:: When using your own data, use any editor (vi(m), gedit, nano, atom, ...) to customize the sample sheet.
.. warning:: **Please ensure not to replace any tabulator symbols with spaces while changing this file.**

We will rewrite this file to fit to the previously downloaded *.fastq.gz* files.



Next, we are going to set up the *config.yaml*.

.. code-block:: bash

    $ cp HRIBO/templates/config.yaml HRIBO/
    $ vi HRIBO/config.yaml

This file contains the following variables:

• **taxonomy** Specify the taxonomic group of the used organism in order to ensure the correct removal of reads mapping to ribosomal genes (Eukarya, Bacteria, Archea). (Option for the preprocessing workflow)
•	**adapter** Specify the adapter sequence to be used. If not set, *Trim galore* will try to determine it automatically. (Option for the preprocessing workflow)
•	**samples** The location of the samples sheet created in the previous step.
•	**genomeindexpath** If the STAR genome index was already precomputed, you can specify the path to the files here, in order to avoid recomputation. (Option for the preprocessing workflow)
•	**uorfannotationpath** If a uORF-annotation file was already pre-computed, you can specify the path to the file here. Please make sure, that the file has the same format as the uORF_annotation_hg38.csv file provided in the git repo (i.e. same number of columns, same column names)
• **alternativestartcodons** Specify a comma separated list of alternative start codons.

.. code-block:: bash

    #Taxonomy of the samples to be processed, possible are Eukarya, Bacteria, Archea
    taxonomy: "Eukarya"
    #Adapter sequence used
    adapter: ""
    samples: "HRIBO/samples.tsv"
    genomeindexpath: ""
    uorfannotationpath: ""
    alternativestartcodons: ""

For this tutorial, we can keep the default values for the *config.yaml*. The organism analyzed in this tutorial is *homo sapiens*, therefore we keep the taxonomy at *Eukarya*. The path to *samples.tsv* is set correctly.

Running the workflow
====================

Now that we have all the required files, we can start running the workflow, either locally or in a cluster environment.

Run the workflow locally
************************

Use the following steps when you plan to execute the workflow on a single server or workstation. Please be aware that some steps
of the workflow might require a lot of memory, specifically for eukaryotic species.

Navigate to the project folder containing the bam/ folder, the annotation.gtf and the genome.fa files and the HRIBO folder. Start the workflow locally from this folder by running:

.. code-block:: bash

    $ snakemake --use-conda -s HRIBO/Snakefile --configfile HRIBO/config.yaml --directory ${PWD} -j 20 --latency-wait 60

Run Snakemake in a cluster environment
**************************************

Use the following steps if you are executing the workflow via a queuing system. Edit the configuration file *cluster.yaml*
according to your queuing system setup and cluster hardware.

Navigate to the project folder on your cluster system. Start the workflow from this folder by running (The following system call shows the usage with Grid Engine.):

.. code-block:: bash

    $ snakemake --use-conda -s HRIBO/Snakefile --configfile HRIBO/config.yaml --directory ${PWD} -j 20 --cluster-config HRIBO/templates/sge-cluster.yaml

Example: Run Snakemake in a cluster environment
***********************************************

.. warning:: **Be advised that this is a specific example, the required options may change depending on your system.**

We ran the tutorial workflow in a cluster environment, specifically a TORQUE cluster environment.
Therefore, we created a bash script *torque.sh* in our project folder.

.. code-block:: bash

    $ vi torque.sh

.. note:: Please note that all arguments enclosed in <> have to be customized. This script will only work if your cluster uses the TORQUE queuing system.
We proceeded by writing the queuing script:

.. code-block:: bash

    #!/bin/bash
    #PBS -N <ProjectName>
    #PBS -S /bin/bash
    #PBS -q "long"
    #PBS -d <PATH/ProjectFolder>
    #PBS -l nodes=1:ppn=1
    #PBS -o <PATH/ProjectFolder>
    #PBS -j oe
    cd <PATH/ProjectFolder>
    source activate HRIBO
    snakemake --latency-wait 600 --use-conda -s HRIBO/Snakefile --configfile HRIBO/config.yaml --directory ${PWD} -j 20 --cluster-config HRIBO/templates/torque-cluster.yaml --cluster "qsub -N {cluster.jobname} -S /bin/bash -q {cluster.qname} -d <PATH/ProjectFolder> -l {cluster.resources} -o {cluster.logoutputdir} -j oe"

We then simply submitted this job to the cluster:

.. code-block:: bash

    $ qsub torque.sh

Using any of the presented methods, this will run the workflow on our dataset and create the desired output files.

Report
******

Once the workflow has finished, we can request an automatically generated *report.html* file using the following command:

.. code-block:: bash

    $ snakemake --use-conda -s HRIBO/Snakefile --configfile HRIBO/config.yaml --report report.html

References
==========

.. bibliography:: references.bib
