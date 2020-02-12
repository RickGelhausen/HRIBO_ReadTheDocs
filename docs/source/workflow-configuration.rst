.. _workflow-configuration:

######################
Workflow configuration
######################

This workflow allows different customization to be able to handle different types of input data.
On this page we explain the different options that can be set to easily customize the workflow.

Default workflow
================

In order to explain what customizations are possible, we will first have a look at the default workflow.

Default:

• Single-end fastq files
• Differential expression analysis: on
• DeepRibo predictions: off

For the default workflow, we expect the fastq files to be in single-end format.
Additionally, we activated differential expression by default. Differential expression requires multiple conditions and RIBO and RNA samples.
A possible *sample.tsv* would look as follows:

+-----------+-----------+-----------+-------------------------+
|   method  | condition | replicate | fastqFile               |
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

By default only reparation predictions are used. The reason for this is that DeepRibo has dependencies that are harder to meet.

No differential expression
==========================

If you do not have multiple conditions and differential expression is activated, you will receive an error message.
To deactivate differential expression, you have to edit the *config.yaml* file.

.. code-block:: bash

    adapter: ""
    samples: "HRIBO/samples.tsv"
    alternativestartcodons: "GTG,TTG"
    # Differential expression: on / off
    differentialexpression: "off"
    # Deepribo predictions: on / off
    deepribo: "off"

This will allow you the use of a sample.tsv like:

+-----------+-----------+-----------+-------------------------+
|   method  | condition | replicate | fastqFile               |
+===========+===========+===========+=========================+
| RIBO      |  A        | 1         | fastq/RIBO-A-1.fastq.gz |
+-----------+-----------+-----------+-------------------------+
| RIBO      |  A        | 2         | fastq/RIBO-A-2.fastq.gz |
+-----------+-----------+-----------+-------------------------+
| RNA       |  A        | 1         | fastq/RNA-A-1.fastq.gz  |
+-----------+-----------+-----------+-------------------------+
| RNA       |  A        | 2         | fastq/RNA-A-2.fastq.gz  |
+-----------+-----------+-----------+-------------------------+

Activating Deepribo
===================

Activating DeepRibo predictions will give you a different file with ORF predictions.
By experience, the top DeepRibo results tend to be slightly better than those of reparation.
For archea, where reparation performs very poorly, DeepRibo is a valid option.

.. note:: In order to use DeepRibo, the tool *singularity* is required. Please refer to the :ref:`overview <overview:Tools>` for details on the installation.

Once you have installed *singularity* turn on DeepRibo in the *config.yaml*:

.. code-block:: bash

    adapter: ""
    samples: "HRIBO/samples.tsv"
    alternativestartcodons: "GTG,TTG"
    # Differential expression: on / off
    differentialexpression: "on"
    # Deepribo predictions: on / off
    deepribo: "on"

When calling snakemake, you will now require additional commandline arguments:

• **--use-singularity:** specify that snakemake can now download and use docker container via singularity.
• **--singularity-args " -c ": specify the *--contain* option to ensure that only the docker containers file system will be used.

If you run deepribo locally
***************************

When running the workflow with DeepRibo locally it might be advised to additionally use the *--greediness 0* option, if you do not have a lot of cores available locally.
This will cause the workflow to submit fewer jobs at the same time. This especially important for DeepRibo as we observed that a single DeepRibo job can finish in less than an hour if it does not have to fight for cores with another DeepRibo job. Otherwise, it can run for several hours at a time.

.. code-block:: bash

    snakemake --use-conda --use-singularity --singularity-args " -c " -s HRIBO/Snakefile --configfile HRIBO/config.yaml --directory ${PWD} -j 10 --latency-wait 60

If you run deepribo on a cluster system
***************************************

When running the workflow with DeepRibo on a cluster system. You have to add the above commandline arguments to your submission script.

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
    snakemake --latency-wait 600 --use-conda --use-singularity --singularity-args " -c " -s HRIBO/Snakefile --configfile HRIBO/config.yaml --directory ${PWD} -j 20 --cluster-config HRIBO/templates/torque-cluster.yaml --cluster "qsub -N {cluster.jobname} -S /bin/bash -q {cluster.qname} -d <PATH/ProjectFolder> -l {cluster.resources} -o {cluster.logoutputdir} -j oe"


.. note:: If you cannot install *singularity* on your cluster, check whether there are modules available for you cluster system.

You can then create an additional submission script that will tell snakemake to activate the module before running jobs.
An example of this would look as follows:

jobscript.sh

.. code-block:: bash

    #!/bin/bash
    module load devel/singularity/3.4.2
    # properties = {properties}
    {exec_job}

Then add the jobscript to the snakemake call:

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
    snakemake --latency-wait 600 --use-conda --use-singularity --singularity-args " -c " --jobscript jobscript.sh -s HRIBO/Snakefile --configfile HRIBO/config.yaml --directory ${PWD} -j 20 --cluster-config HRIBO/templates/torque-cluster.yaml --cluster "qsub -N {cluster.jobname} -S /bin/bash -q {cluster.qname} -d <PATH/ProjectFolder> -l {cluster.resources} -o {cluster.logoutputdir} -j oe"

This will specify to snakemake that it will execute *module load devel/singularity/3.4.2* when submitting each job.

.. note:: This is a specific example for our TORQUE cluster system. The specific way of loading modules, as well as the available modules, can differ on each system.


Paired-end support
==================

We allow paired-end data in our workflow.
Unfortunately, many of the downstream tools, like the prediction tools, cannot use paired-end data.
Therefore, we use the tool *flash2* **TODO cite/link** to convert paired-end data to single-end data.

In order to use paired-end data, simply replace the *Snakefile* with the *Snakefile_pairedend*.
This will now require a special *samples_pairedend.tsv*, which is also available in the HRIBO templates folder.

+-----------+-----------+-----------+----------------------------+----------------------------+
|   method  | condition | replicate | fastqFile1                 | fastqFile2                 |
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
