.. _faq:

##########################
Frequently asked questions
##########################

Q: When using singularity I get ``ERROR  : Failed to set effective UID to 0``.
==============================================================================

Prior to the installation of singularity change ``with_suid=1`` to ``with_suid=0`` in the mconfig file in the singularity folder.
This should not be necessary for newer versions of singularity.
