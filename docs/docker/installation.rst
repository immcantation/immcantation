Getting the /suite Container
--------------------------------------------------------------------------------

Requires an installation of Docker 1.9+ or Singularity 2.3+.

Docker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. parsed-literal::

    # Pull release version |docker-version|
    docker pull immcantation/suite:|docker-version|

    # Pull the latest development build
    docker pull immcantation/suite:devel


Our containers are Linux-based, so if you are using a Windows computer,
please make sure that you are using Linux containers and not Windows containers
(this can be changed in `Docker Desktop <https://www.docker.com/products/docker-desktop/>`_
and won't affect your existing containers).


Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. parsed-literal::

    # Pull release version |docker-version|
    IMAGE="immcantation_suite-|docker-version|.sif"
    singularity build $IMAGE docker://immcantation/suite:|docker-version|

The instructions to use containers from `Docker Hub <https://hub.docker.com/>`_
with Singularity can be slightly different for different versions of Singularity.
If the command shown above doesn't work for you, please visit
`Singularity Documentation <https://www.sylabs.io/docs/>`_ and look for the
specific command for your Singularity version under *Build a container*.
