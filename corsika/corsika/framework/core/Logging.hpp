/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/**
 *  @file Logging.hpp
 *
 * CORSIKA8 logging utilities.
 *
 * See testLogging.cpp for a complete set of examples for
 * how the logging functions should be used.
 */

#pragma once

// Configure some behaviour of sdlog.
// This must be done before spdlog is included.

// use the coarse system clock. This is *much* faster
// but introduces a timestamp error of O(10 ms) which is fine for us.
#ifndef SPDLOG_CLOCK_COARSE
#define SPDLOG_CLOCK_COARSE
#endif
// do not create a default logger (we provide our own "corsika" logger)
#ifndef SPDLOG_DISABLE_DEFAULT_LOGGER
#define SPDLOG_DISABLE_DEFAULT_LOGGER
#endif
// use __PRETTY_FUNCTION__ instead of __FUNCTION__ where
// printing function names in trace statements. This is much
// nicer than __FUNCTION__ under GCC/clang.
#ifndef SPDLOG_FUNCTION
#define SPDLOG_FUNCTION __PRETTY_FUNCTION__
#endif

// if this is a Debug build, include debug messages in objects
#ifdef _C8_DEBUG_
// trace is the highest level of logging (ALL messages will be printed)
#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif
#else // otherwise, remove everything but "error" and worse messages
#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif
#endif

#include <spdlog/fmt/ostr.h> // will output whenerver a streaming operator is found
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace corsika {

  /*
   * The default pattern for CORSIKA8 loggers.
   */
  const std::string minimal_pattern{"[%n:%^%-8l%$] %v"};
  const std::string default_pattern{"[%n:%^%-8l%$(%s:%#)] %v"};
  const std::string source_pattern{"[%n:%^%-8l%$(%s:%!:%#)] %v"};

  /**
   * Create a new C8-style logger.
   *
   * Use this if you are explicitly (and can guarantee) that you
   * are creating a logger for the first time. It is recommended
   * that for regular usage, the `get_logger` function is used instead
   * as that will also create the logger if it has not yet been created.
   *
   * Calling `create_logger` twice to create the same logger will
   * result in an spdlog duplicate exception.
   *
   * @param name           The unique name of the logger.
   * @param defaultlog     If True, set this as the default logger.
   * @returns              The constructed and formatted logger.
   */
  std::shared_ptr<spdlog::logger> create_logger(std::string const& name,
                                                bool const defaultlog = false);

  /**
   * Get a smart pointer to an existing logger.
   *
   * This should be the default method for code to obtain a
   * logger. If the logger *does not* exist, it is *created* and
   * returned to the caller.
   *
   * This should be preferred over `create_logger`.
   *
   * @param name    The name of the logger to get.
   * @param defaultlog   If True, make this the default logger.
   * @returns              The constructed and formatted logger.
   */
  std::shared_ptr<spdlog::logger> get_logger(std::string const& name,
                                             bool const defaultlog = false);

  /**
   * The default "corsika" logger.
   */
  static inline std::shared_ptr<spdlog::logger> corsika_logger =
      get_logger("corsika", true);

  // many of these free functions are special to the logging
  // infrastructure so we hide them in the corsika::logging namespace.
  namespace logging {

    // bring spdlog into the corsika::logging namespace
    using namespace spdlog;

    /**
     * Set the default log level for all *newly* created loggers.
     *
     *  @param minlevel The minimum log level required to print.
     *
     */
    auto set_default_level(level::level_enum const minlevel) -> void;

    /**
     * Add the source (filename, line no) info to the logger.
     *
     * @param logger  The logger to set the level of.
     *
     */
    template <typename TLogger>
    auto add_source_info(TLogger& logger) -> void;

    /**
     * Reset the logging pattern to the default.
     *
     * @param logger  The logger to set the level of.
     *
     */
    template <typename TLogger>
    auto reset_pattern(TLogger& logger) -> void;

  } // namespace logging

// define our macro-style loggers
// these use the default "corsika" logger
#define CORSIKA_LOG_TRACE SPDLOG_TRACE
#define CORSIKA_LOG_DEBUG SPDLOG_DEBUG
#define CORSIKA_LOG_INFO SPDLOG_INFO
#define CORSIKA_LOG_WARN SPDLOG_WARN
#define CORSIKA_LOG_ERROR SPDLOG_ERROR
#define CORSIKA_LOG_CRITICAL SPDLOG_CRITICAL

// and the specific logger versions
// these take a logger instance as their first argument
#define CORSIKA_LOGGER_TRACE SPDLOG_LOGGER_TRACE
#define CORSIKA_LOGGER_DEBUG SPDLOG_LOGGER_DEBUG
#define CORSIKA_LOGGER_INFO SPDLOG_LOGGER_INFO
#define CORSIKA_LOGGER_WARN SPDLOG_LOGGER_WARN
#define CORSIKA_LOGGER_ERROR SPDLOG_LOGGER_ERROR
#define CORSIKA_LOGGER_CRITICAL SPDLOG_LOGGER_CRITICAL

} // namespace corsika

#include <corsika/detail/framework/core/Logging.inl>
#include <corsika/detail/framework/core/SpdlogSpecializations.inl>
