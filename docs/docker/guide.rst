Using the Container
================================================================================

Invoking a shell inside the container
--------------------------------------------------------------------------------

To invoke a shell session inside the container:

.. parsed-literal::

    # Docker command
    docker run -it kleinstein/immcantation:|docker-version| bash

    # Singularity command
    singularity shell immcantation-|docker-version|.img

Sharing files with the container
--------------------------------------------------------------------------------

Sharing files between the host operating system and the container requires you
to bind one of the container's mount points to a folder on the host using the
``-v`` argument to ``docker`` or the ``-B`` argument to ``singularity``.
There are four available mount points defined in the container::

    /data
    /scratch
    /software
    /oasis

To invoke a shell session inside the container with ``$HOME/project`` mounted to
``/data``:

.. parsed-literal::

    # Docker command
    docker run -it -v $HOME/project:/data:z kleinstein/immcantation:|docker-version| bash

    # Singularity command
    singularity shell -B $HOME/project:/data immcantation-|docker-version|.img

Note, the ``:z`` in the ``-v`` argument of the ``docker`` command is essential.

Executing a specific command
--------------------------------------------------------------------------------

To execute a specific command inside the container with ``$HOME/project`` mounted to
``/data``:

.. parsed-literal::

    # Docker command
    docker run -v $HOME/project:/data:z kleinstein/immcantation:|docker-version| versions report

    # Singularity command
    singularity exec -B $HOME/project:/data immcantation-|docker-version|.img versions report

In this case, we are executing the ``versions report`` command which will inspect
the installed software versions and print them to standard output.

There is an analagous ``builds report`` command to display the build date and
changesets used during the image build. This is particularly relevant if you
are using the ``kleinstein/immcantation:devel`` development builds.
