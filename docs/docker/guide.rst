.. _DockerGuide:

Using the Container
================================================================================

Invoking a shell inside the container
--------------------------------------------------------------------------------

To invoke a shell session inside the container:

.. parsed-literal::

    # Docker command
    docker run -it immcantation/suite:|docker-version| bash

    # Singularity command
    singularity shell immcantation_suite_|docker-version|.sif

Sharing files with the container
--------------------------------------------------------------------------------

Sharing files between the host operating system and the container requires you
to bind a directory on the host to one of the container's mount points using the
``-v`` argument for ``docker`` or the ``-B`` argument for ``singularity``.
There are four available mount points defined in the container::

    /data
    /scratch
    /software
    /oasis

For example, to invoke a shell session inside the container with ``$HOME/project`` mounted to
``/data``:

.. parsed-literal::

    # Docker command
    docker run -it -v $HOME/project:/data:z immcantation/suite:|docker-version| bash

    # Singularity command
    singularity shell -B $HOME/project:/data immcantation_suite_|docker-version|.sif

Note, the ``:z`` in the ``-v`` argument of the ``docker`` command is essential.

Executing a specific command
--------------------------------------------------------------------------------

After invoking an interactive session inside the container, commands can be
executed in the container shell as they would be executed in the host shell.

Alternatively, it is possible to execute a specific command directly inside the 
container without starting an interactive session. The next example demonstrates 
how to execute ``versions report`` with ``$HOME/project`` mounted to ``/data``:

.. parsed-literal::

    # Docker command
    docker run -v $HOME/project:/data:z immcantation/suite:|docker-version| versions report

    # Singularity command
    singularity exec -B $HOME/project:/data immcantation_suite_|docker-version|.sif versions report

In this case, we are executing the ``versions report`` command which will inspect
the installed software versions and print them to standard output.

There is an analagous ``builds report`` command to display the build date and
changesets used during the image build. This is particularly relevant if you
are using the ``immcantation/suite:devel`` development builds.
