Contributing to CORSIKA 8
=========================

The CORSIKA Project welcomes contributions. This guide explains how to find the
right place to contribute and how work is integrated into the project. Connect
to `GitLab <https://gitlab.iap.kit.edu>`_ and the
`corsika-devel mailing list <https://www.lists.kit.edu/sympa/subscribe/corsika-devel>`_
to get in touch with the project.

Contribution Workflow
---------------------

* Development is organized through GitLab issues for feature requests, ideas,
	discussions, and bug fixes.
* Issues are discussed within the project, assigned to milestones, and worked
	on in branches.
* Proposed changes are reviewed and merged through merge requests.
* Every merge request undergoes code review and must be approved before merge.
	The project wiki defines the detailed
	`code-approval procedure <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/Code-Approval-Procedure>`_.
* For branch, release, and versioning policy, see :doc:`gitflow`.
* For additional branch-update background, see :doc:`git_workflow`.

Code Formatting and Conventions
-------------------------------

CORSIKA 8 uses ``clang-format`` to keep source formatting consistent. Integrate
it with your editor when possible, and run the project's formatting script
before committing. The source-format definition is maintained in
``.clang-format``. See :doc:`coding_conventions` and :doc:`unit_tests` for the
detailed project rules.

Use ``do-clang-tidy.py`` for static analysis of the C++ coding conventions. It
requires a CMake build directory containing ``compile_commands.json``; configure
CMake with ``-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`` when necessary. The default
checks translation units changed relative to ``main``:

.. code-block:: console

	$ python3 do-clang-tidy.py --build-dir build

Use ``--all`` in CI to check every configured translation unit. The script
treats diagnostics as errors. Fix the code or use a narrowly scoped ``NOLINT``
comment with a justification when a specific finding is a deliberate exception.
The project policy is maintained in ``.clang-tidy``.

For new physics functionality, start with
:doc:`../framework/howto_create_module` after reading these project
requirements.

AI-assisted contributions must also follow :doc:`ai_use`.

Release Versioning
------------------

The full release workflow and version-number policy are defined in
:doc:`gitflow`.

Scientific Authorship
---------------------

The CORSIKA Project decides who becomes a scientific author. Relevant criteria
include responsibility for a functionality or management area, active
participation in the project, adherence to this guide, agreement to the
`usage and collaboration agreement <https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/blob/main/usage_and_collaboration_agreement.md>`_,
and agreement among project members.