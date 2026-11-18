Unit Tests
==========

CORSIKA unit tests use the
`Catch2 framework <https://github.com/catchorg/Catch2/blob/devel/docs/assertions.md>`_.
Prefer the ``CHECK`` family of assertions because it continues a test after a
failure; use ``REQUIRE`` only when continuing would make the test invalid.

Test Scope
----------

Test every public interface and its basic behavior. Each functional group of
classes should have an associated unit test. New functionality must include
tests, and merge requests must not reduce the measured code coverage. The
project aims for coverage to approach 100 percent over time.

Tests must be fast, robust, and independent:

* Avoid extensive physics simulations and statistical tests with many events.
* Check expected behavior and formulas with the smallest useful input.
* Replace complicated input or output objects with simplified test doubles when
  that isolates the behavior under test.
* Make each test depend on as little code and shared state as possible.

These expectations are enforced together with the project CI requirements
described in :doc:`contributing`.