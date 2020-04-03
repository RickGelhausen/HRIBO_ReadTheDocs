.. _faq:

##########################
Frequently asked questions
##########################

Q: This could avoid potential errors like: ``ERROR  : Failed to set effective UID to 0``.
=========================================================================================

A: Prior to the installation of singularity change ``with_suid=1`` to ``with_suid=0`` in the mconfig file in the singularity folder.
