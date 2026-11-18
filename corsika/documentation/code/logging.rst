Logging
========

CORSIKA provides a logging interface built around
`spdlog <https://github.com/gabime/spdlog>`_, which uses
`fmt <https://fmt.dev/latest/index.html>`_ for typed message formatting.

Logging Interfaces
------------------

Two interfaces are available:

* Runtime logging uses functions such as ``logging::info("message")`` or
  methods such as ``logger->info("message")``.
* Compile-time-removable logging uses macros such as ``CORSIKA_LOG_INFO`` and
  ``CORSIKA_LOGGER_INFO``. These calls can be removed in builds where frequent
  logging would affect performance.

Calls without a specific logger use the global ``corsika`` logger. Reserve that
logger for framework code. Every physics process should create and retain its
own logger with a lower-case process name:

.. code-block:: c++

   auto logger{get_logger("process_name")};
   logger->info("Process initialized with {} particles.", particle_count);
   CORSIKA_LOGGER_DEBUG(logger, "Detailed process state: {}", state);

Use the logger instance rather than the global free functions in process code.
The spdlog `usage samples <https://github.com/gabime/spdlog#usage-samples>`_
contain further formatting examples.

Choosing a Level
----------------

The available levels, ordered by severity, are ``trace``, ``debug``, ``info``,
``warn``, ``error``, and ``critical``.

* Use ``info`` for infrequent, non-essential run information such as
  initialization parameters or summary statistics.
* Use ``warn`` for infrequent parameter states that are outside an appropriate
  range but do not stop the simulation.
* Use ``error`` for recoverable failures, for example when one shower or
  particle must be stopped while the simulation can continue.
* Use ``critical`` immediately before an exception terminates the simulation for
  an unrecoverable failure.
* Use the ``CORSIKA_LOG_DEBUG`` or ``CORSIKA_LOGGER_DEBUG`` macros for detailed
  messages that may occur a few times per event.
* Use the ``CORSIKA_LOG_TRACE`` or ``CORSIKA_LOGGER_TRACE`` macros for
  per-particle or otherwise high-frequency messages.

Dynamic Log Levels
------------------

Set a logger's runtime threshold with ``set_level``:

.. code-block:: c++

   logger->set_level(logging::level::critical);

Choose from ``critical``, ``error``, ``warning``, ``info``, ``debug``, and
``trace``. Compile-time-removable macro calls remain unavailable when the build
configuration removes them, so enabling a runtime level cannot restore those
messages without rebuilding the framework.