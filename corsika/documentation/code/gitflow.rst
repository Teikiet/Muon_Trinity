.. _release-versioning:

Gitflow and Release Versioning
==============================

This document defines the CORSIKA 8 development and release workflow. It uses
Gitflow: ``main`` contains released code, while versioned ``dev/X.Y.0`` branches
are the shared integration branches for planned minor releases. Naming
development branches after their target release allows more than one release
line to remain active when a feature needs a longer development timeline.

Long-lived Branches
-------------------

* ``main`` contains production-ready releases only. Every release commit is
  tagged ``vX.Y.Z``; direct pushes are not permitted.
* ``dev/X.Y.0`` contains the integrated work targeted at release ``X.Y.0``.
  Feature merge requests target the appropriate versioned development branch
  until its feature freeze; direct pushes are not permitted.
* ``release/X.Y.Z`` is the read-only archive of a released development line. It
  is created by renaming ``dev/X.Y.0`` after release ``X.Y.Z`` is tagged on
  ``main``. This should allow better backport of bugfixes.

``main`` and every active ``dev/X.Y.0`` branch are protected and require a
passing pipeline and the review requirements in :doc:`contributing`. Every
``release/X.Y.Z`` branch is protected in read-only mode.

Short-lived Branches
--------------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 45

   * - Branch
     - Created from
     - Merges into
     - Purpose
   * - ``feature/<issue>-<slug>``
     - ``dev/X.Y.0``
     - ``dev/X.Y.0``
     - One issue, feature, refactoring, or ordinary bug fix.
   * - ``release/X.Y.Z``
     - renamed from ``dev/X.Y.0``
     - none
     - Read-only archive of a released version.
   * - ``hotfix/vX.Y.Z``
     - ``main``
     - ``main``, then active ``dev/X.Y.0`` branches
     - Repair an already released version.

Use lowercase kebab-case for ``<slug>``. Keep a branch scoped to one issue and
delete it after its merge request is accepted. Do not merge a feature branch
directly into ``main``.

.. code-block:: text

   gitGraph
     commit id: "release v8.0.0" tag: "v8.0.0"
     branch dev/8.1.0
     checkout dev/8.1.0
     commit id: "integrate"
     branch feature/123-example
     checkout feature/123-example
     commit id: "implement"
     checkout dev/8.1.0
     merge feature/123-example id: "reviewed feature"
     commit id: "physics validation"
     checkout main
     merge dev/8.1.0 tag: "v8.1.0"

After the tagged release, rename ``dev/8.1.0`` to ``release/8.1.0`` and protect
it as read-only. The Mermaid ``gitGraph`` notation cannot represent a branch
rename; the final rename is therefore described in the release flow below.

Development Flow
----------------

#. Discuss and assign the GitLab issue, then select the target release and
   create ``feature/<issue>-<slug>`` from its ``dev/X.Y.0`` branch.
#. Implement one cohesive change, add tests and documentation, and keep the
   branch current by merging its target ``dev/X.Y.0`` branch when needed.
#. Open a merge request from the feature branch to the target ``dev/X.Y.0``
   branch. It must pass CI and complete code review before merge.
#. Delete the feature branch after the merge request is accepted.

Rebasing and Merging
--------------------

Merge the current target branch into a feature branch whenever it needs to pick
up integration changes. This is the default for long-lived or shared feature
branches because it preserves published history and records the integration
point explicitly:

.. code-block:: shell

   git fetch origin
   git switch feature/123-example
   git merge origin/dev/8.2.0
   # Resolve conflicts, run the relevant tests, then create the merge commit.
   git push

Rebasing is optional for a short-lived feature branch owned by one developer,
typically before opening its merge request when a linear history makes review
easier. It rewrites commit identities, so never rebase a shared branch. After
rebasing a published personal branch, use ``git push --force-with-lease``, never
plain ``--force``.

Never rebase, force-push, or push directly to ``main``, ``dev/X.Y.0``, or
``release/X.Y.Z``. These protected branches are updated only through approved
merge requests. The GitLab merge strategy configured for the project determines
whether an accepted feature merge request produces a merge commit or a squash
merge; contributors must not merge locally into a protected branch.

Use ordinary merges, not rebases, when propagating a released or hotfixed
change from ``main`` into active ``dev/X.Y.0`` branches. These merges preserve
the shared history and make it clear which release changes have reached each
development line.

Physics Validation and Feature Freeze
-------------------------------------

For each planned release, the steering board selects and announces a feature
freeze date for the corresponding ``dev/X.Y.0`` branch. Until that date, feature
merge requests may target the branch. From the freeze onward, only bug fixes
and their necessary regression tests and documentation changes may be merged.
New features, physics models, interface changes, and unrelated refactoring are
deferred to a later ``dev/X.Y.0`` branch.

The steering board also agrees the timeframe in which each release candidate
must become physics complete. Physics validation is performed on the frozen
development branch and establishes that the agreed physics scope is complete
and validated. A feature or part of a feature that is not physics complete by
the agreed deadline must not be included in the release.

* If the work has not merged, retarget its feature branch and merge request to
  the next appropriate ``dev/X.Y.0`` branch.
* If the work has already merged into the frozen branch, create a reviewed
  revert merge request for ``dev/X.Y.0`` so the release line remains physics
  complete. Then cherry-pick the feature's reviewed commits onto the next
  appropriate ``dev/X.Y.0`` branch through a separate merge request, resolving
  conflicts and rerunning validation there.

Do not merge the frozen development branch into the later development branch to
defer an incomplete feature: that would also carry release-only fixes and other
changes not intended for the later release line. Record the original merge
request and the revert and cherry-pick merge requests in the issue so the
feature's history remains traceable.

After the steering board accepts the completed physics validation, merge the
frozen ``dev/X.Y.0`` branch into ``main``, create the annotated release tag
``vX.Y.Z``, then rename the development branch to ``release/X.Y.Z`` and protect
it in read-only mode. The archived branch preserves the exact development line
used for the release; new work belongs on an active ``dev/X.Y.0`` branch.

For a released-version defect, create ``hotfix/vX.Y.Z`` from the affected
``main`` release. Validate it as a release, merge it into ``main``, tag it, and
merge the resulting ``main`` commit into every active ``dev/X.Y.0`` branch.
Hotfixes must never wait for unrelated work already in a development branch.

Create the next ``dev/X.(Y+1).0`` branch from the current development line when
work intended for the following minor release must proceed before ``X.Y.0`` is
released. A feature branch targets the earliest versioned development branch in
which it can safely be released. This allows a long-running feature to continue
on a later ``dev/X.Y.0`` branch while an earlier release line is stabilized.

Version Numbers
---------------

CORSIKA 8 releases use ``X.Y.Z``, as recorded in the top-level CMake project
version. This is a project-specific compatibility policy rather than a claim of
Semantic Versioning:

* ``X`` remains ``8`` for the foreseeable future. Increment it only for a
  massive overhaul that establishes a new CORSIKA baseline.
* Increment ``Y`` for approved features or bug fixes that may include interface
  changes.
* Increment ``Z`` for compatible bug fixes or features that do not change
  public interfaces. Production code remains interface-compatible across ``Z``
  updates.

Set the planned release version in ``CMakeLists.txt`` during preparation of the
corresponding ``dev/X.Y.0`` branch, or on a ``hotfix/`` branch. Do not increment
the version speculatively on feature branches. Final releases receive the
annotated Git tag ``vX.Y.Z`` on ``main``.

When a release candidate is needed, use a pre-release tag such as
``vX.Y.Z-rc.1``; the CMake project version remains the numeric ``X.Y.Z`` until
the final release.

Merge Request Expectations
--------------------------

All merge requests follow :doc:`contributing`. In addition:

* Feature merge requests target their selected ``dev/X.Y.0`` branch and link
  their issue.
* The merge request from a frozen ``dev/X.Y.0`` branch to ``main`` identifies
  the target version and includes the version update, release notes, physics
  validation evidence, test evidence, and required documentation.
* Hotfix merge requests identify the target version and include the version
  update, test evidence, and required documentation.
* ``main`` accepts merges only from frozen ``dev/X.Y.0`` branches and
  ``hotfix/`` branches.
* Every merge from ``main`` is followed by a merge request or approved merge
  back into each active ``dev/X.Y.0`` branch that needs the change.
