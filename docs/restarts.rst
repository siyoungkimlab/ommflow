Restarts and wall time
======================

Automatic restarts
------------------

A missing or empty ``workdir`` starts a new run, so a directory pre-created by
a batch scheduler is not mistaken for an unresumable one. An existing
``workdir`` holding a previous run is automatically treated as a restart and
must contain ``solvated.pdb``, ``system.xml``, ``integrator.xml``,
``checkpoint.chk``, and ``final.toml``.

For a wall-time-limited job, submit the same command again:

.. code-block:: console

   ommflow protein.pdb \
     --workdir protein_md \
     --production-ns 200

The saved checkpoint restores production state and the workflow appends to
``trajectory.dcd`` and ``state.csv``.

Production target
-----------------

``production_ns`` is an absolute production-time target, not a per-submission
duration. Each submission runs only the remaining time required to reach the
target. If ``--production-ns`` is omitted during a restart, the target saved
in ``workdir/final.toml`` is used.

Extend a completed or ongoing target with an absolute value:

.. code-block:: console

   ommflow \
     --workdir protein_md \
     --production-ns 300

This updates ``final.toml`` and runs only the remaining production time.

Restored settings
-----------------

A restart reloads every setting saved in ``workdir/final.toml`` — the
equilibration, reporting, and checkpoint intervals, the platform, and all
early-stop settings — unless the same setting is supplied again on the command
line or in ``--config``. This keeps ``trajectory.dcd`` and ``state.csv`` on one
cadence across submissions.

``final.toml`` is rewritten only after a restart passes every check, so a
rejected resume leaves the run directory exactly as it was and still usable.

Reporting and checkpoints
-------------------------

Equilibration and production have separate reporting intervals. For
restart-safe trajectory appending, ``checkpoint_interval_ns`` must be no
greater than, and divide, ``production_report_interval_ns`` exactly in
integration steps:

.. code-block:: toml

   production_report_interval_ns = 1.0
   checkpoint_interval_ns = 0.01

Use the same OpenMM version and, where possible, the same execution platform
when restoring a checkpoint.

Early-stop restarts
-------------------

For an early-stop run, ``pocket.json`` is created after equilibration and is
never recomputed or re-resolved on restart. Its saved target and production
atom mapping remain authoritative. ``monitor.csv`` appends without another header,
and ``status.json`` retains the consecutive-detachment counter. A confirmed
``detached`` status is terminal when rerunning with the same target and
``early_stop`` enabled. To deliberately continue, disable early stopping or
provide a larger explicit ``--production-ns`` target. A ``target_reached``
status safely does nothing until the absolute target is increased.

Do not enable early stopping for an already-existing work directory that lacks
``pocket.json``: start a new work directory so the production-start pocket can
be defined correctly.
