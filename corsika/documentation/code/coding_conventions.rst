Coding Conventions
==================

These conventions keep the CORSIKA 8 codebase understandable, maintainable,
and consistent. They supplement the project's ``clang-format`` configuration
and apply to framework code and modules alike.

Project Structure
-----------------

The source tree is organized relative to the CMake project source directory:

* ``corsika`` contains public, mostly header-only framework code.
* ``modules`` contains additional physics packages integrated with the framework.
* ``src`` contains framework code that is compiled into libraries.
* ``examples`` contains executable examples and ``tests`` contains unit tests.
* ``cmake``, ``externals``, ``documentation``, and ``tools`` contain build,
  dependency, documentation, and supporting resources.

All framework facilities are available in the ``corsika`` namespace. Modules may
use their own nested namespaces where appropriate. Code should be sufficiently
clear for a new contributor to understand and extend.

C++ Files and Naming
--------------------

* Use ``.hpp`` for public headers, ``.cpp`` for executable or library entry
  points, and ``.inl`` for inline implementations included by the corresponding
  header.
* Headers declare and document interfaces. Keep implementation details in
  ``.inl`` files, and do not place Doxygen documentation in ``.inl`` files.
* Guard files intended for inclusion with ``#pragma once``. Include project
  headers as ``<path/to/file.hpp>``.
* A class and its containing file use the same PascalCase name, such as
  ``class Doer`` in ``Doer.hpp``.
* Use ``snake_case`` for free functions and namespaces, PascalCase for classes
  and enum types, camelCase for methods, and camelCase with a trailing
  underscore for data members.
* Avoid abbreviations where a descriptive name is practical. In particular, do
  not use single-letter loop variables beyond a single short expression.

Declarations and Class Design
-----------------------------

Place ``const`` at the right of the type, for example ``Type const& value``.
Apply the same principle to ``constexpr`` declarations. Define enums as locally
and strongly scoped as practical, preferably in their owning class.

Classes should expose behavior rather than public data. Keep implementation
details private or protected, use ``const`` for operations that do not modify
state, and provide accessors only where they add a useful interface. Prefer
positive state predicates such as ``isValid()`` or ``hasData()``. Declare data
members after methods.

Consider whether a class owns resources and whether it should be movable or
copyable. Follow the rule of zero or five: use compiler defaults where possible;
when resource management requires custom special members, define the complete
set of copy and move operations deliberately. Use ``= default`` or ``= delete``
to state the intended behavior explicitly.

Use C++23 facilities. Keep type aliases near the top of the class or method that
uses them, name aliases in lower ``snake_case`` with a ``_type`` suffix, and
prefix template type parameters with ``T``. Provide deduction guides and
free-function instance makers for class templates when they make type deduction
clearer.

C++23 Facilities
----------------

Prefer language and library facilities that make interfaces explicit without
introducing unnecessary abstractions. In particular:

* Use concepts and ``requires`` clauses to state template requirements instead
  of SFINAE-based enablement or hand-written detection traits. Keep concepts
  small and name them for the capability they describe.
* Use ``std::span`` for a non-owning view of contiguous input. Do not use it to
  express ownership or retain it after the referenced storage can disappear.
* Prefer ``std::ranges`` algorithms and views when they make a transformation
  or selection clearer than an iterator-based equivalent.
* Use ``std::expected`` for an anticipated, recoverable result that callers are
  expected to handle. Continue to use exceptions for rare, major failures.
* Mark results that must be checked with ``[[nodiscard]]``.

General Rules
-------------

* Do not use ``std::cout``, ``std::cerr``, or ``std::endl``. Use the project
  :doc:`logging` interface instead.
* Build with warnings enabled. Fix warnings originating in CORSIKA code; if a
  third-party warning must be suppressed locally, document the reason at the
  suppression site with ``FIXME``.
* All unit tests must pass on the configured CI systems. New code must include
  appropriate tests and must not reduce coverage. See :doc:`unit_tests`.
* Throw exceptions for rare, major failures. Do not use exception catching as a
  routine error-handling mechanism without an agreed special case.
* Prefer references where appropriate. Do not expose raw pointers in public
  interfaces; prefer standard smart pointers when ownership is required.
* New physics algorithms require appropriate validation in addition to unit
  tests.

Simple Examples
---------------

The following small examples illustrate the conventions above. They are not
complete physics implementations.

An ``EnergySelector.hpp`` header can use a PascalCase class name, a snake_case
namespace, camelCase methods, and a trailing underscore for its data member:

.. code-block:: c++

   #pragma once

   namespace corsika::example_module {

   class EnergySelector {
   public:
     explicit EnergySelector(double const minimumEnergy)
         : minimumEnergy_(minimumEnergy) {}

     [[nodiscard]] bool accepts(double const energy) const {
       return energy >= minimumEnergy_;
     }

   private:
     double minimumEnergy_;
   };

   } // namespace corsika::example_module

Document the public interface by describing its purpose, parameters, and return
value instead of repeating the implementation:

.. code-block:: c++

   /**
    * Checks whether an energy passes the configured threshold.
    *
    * \param energy Energy to check.
    * \return ``true`` when ``energy`` is at least the threshold.
    */
   [[nodiscard]] bool accepts(double const energy) const;

Keep the corresponding CMake declaration small, consistently formatted, and
explicit about the target it creates:

.. code-block:: cmake

   add_library (example_module
     EnergySelector.cpp
   )

   target_link_libraries (example_module
     PUBLIC corsika
   )

Comments and Doxygen
--------------------

Public C++ declarations must have useful Doxygen documentation in their header
files. Comments should clarify non-obvious intent, constraints, or trade-offs;
remove commented-out debugging code and avoid restating code that is already
self-explanatory.

Start a Doxygen block with ``/**``, ``///!``, or ``/*!``. The first paragraph is
the brief description, so ``\brief`` is normally unnecessary. Doxygen comments
apply to the following declaration, so do not add redundant ``\class`` or
``\function`` commands. Use tags such as ``\param``, ``\return``, ``\sa``,
``\ref``, ``\example``, ``\todo``, and ``\bug`` only when they add useful
information. Add a ``\file`` block only when a file itself documents global
variables, enums, or free functions.

Other Languages and CMake
-------------------------

Use ``.f`` for Fortran, ``.py`` for Python, and ``.h`` and ``.c`` for C headers
and source files. Follow established local style for CMake:

* Write CMake commands in lower case and variables in upper case.
* Put a space before an opening parenthesis, for example ``set (VARIABLE value)``.
* Use two spaces for logical indentation, wrap long invocations by keyword, and
  choose expressive variable and function names.
* Comment CMake logic where its purpose is not obvious from the code.