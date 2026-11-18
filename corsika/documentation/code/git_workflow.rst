Git Workflow: Rebase and Merge
==============================

Feature branches can drift away from ``main`` while development continues. To
test a feature against the latest project state, integrate the current upstream
branch into your feature branch deliberately.

For the authoritative CORSIKA 8 branch, release, and versioning policy, see
:doc:`gitflow`.

Both ``git rebase`` and ``git merge`` are useful, but they serve different
purposes.

.. figure:: https://gitlab.iap.kit.edu/AirShowerPhysics/corsika/-/wikis/uploads/0a330dc12d8a281abf21464626d2c2df/image.png
   :alt: Schematic comparison of git merge and git rebase histories

   Schematic comparison of merge and rebase histories.

Rebase
------

Rebasing creates a flat and clean history by replaying your feature-branch
commits on top of the current upstream branch. This is often useful when you
work alone on a branch and want to test or continue the feature against the
latest project state.

The important caveat is that rebasing rewrites commit history. Commits after the
rebase are different commits from the ones before the rebase. This can create
serious problems when several people have checked out the same branch.

Merge
-----

Merging preserves the actual development history and creates a merge commit.
This is appropriate when several people actively work on the same branch, when
coordination around a rewritten branch is impractical, or when merging a feature
branch back into the main project history.

Frequent merges from ``main`` into a long-lived feature branch can make the
feature branch noisy, so avoid long-lived broad feature branches where possible.

CORSIKA 8 Guidelines
--------------------

* Keep feature branches focused and short-lived where possible.
* If you work alone, use ``git rebase main`` or ``git rebase master`` to place
  your branch on top of the current upstream branch before final testing.
* If other people have checkouts of your branch, inform them before rebasing and
  coordinate how they will update to the rewritten history.
* If such coordination is not practical, use ``git merge`` instead of rebasing.
* When closing a feature branch and integrating it into the main project, use a
  merge or fast-forward merge as suggested by GitLab.
* Squashing can be useful in special cases, but decide this during review. It
  loses detailed commit history, which is not always desirable.

Further Reading
---------------

* `Why you should stop using git rebase <https://medium.com/@fredrikmorken/why-you-should-stop-using-git-rebase-5552bee4fed1>`_
* `Git rebase vs merge <https://medium.com/datadriveninvestor/git-rebase-vs-merge-cc5199edd77c>`_
