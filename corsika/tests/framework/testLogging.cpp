/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#include <corsika/framework/core/Logging.hpp>

#include <catch2/catch_all.hpp>

using namespace corsika;

TEST_CASE("Logging", "[Logging]") {

  SECTION("top level functions using default corsika logger") {
    logging::info("(1) This is an info message!");
    logging::warn("(1) This is a warning message!");
    logging::debug("(1) This is a debug message!");
    logging::error("(1) This is an error message!");
    logging::critical("(1) This is a critical error message!");
  }

  SECTION("create a specific logger") {

    // create a logger manually
    auto logger = create_logger("loggerA");

    // set a custom pattern for this logger
    logger->set_pattern("[%n:%^%-8l%$] custom pattern: %v");

    // and make sure we can log with this created object
    logger->info("(2) This is an info message!");
    logger->warn("(2) This is a warning message!");
    logger->debug("(2) This is a debug message!");
    logger->error("(2) This is an error message!");
    logger->critical("(2) This is a critical error message!");

    // get a reference to the logger using Get
    auto other = get_logger("loggerA");

    // and make sure we can use this other reference to log
    other->info("(3) This is an info message!");
    other->warn("(3) This is a warning message!");
    other->debug("(3) This is a debug message!");
    other->error("(3) This is an error message!");
    other->critical("(3) This is a critical error message!");
  }

  SECTION("get a new logger") {

    // get a reference to an unknown logger
    auto logger = get_logger("loggerB");

    // and make sure we can log with this created object
    logger->info("(4) This is an info message!");
    logger->warn("(4) This is a warning message!");
    logger->debug("(4) This is a debug message!");
    logger->error("(4) This is an error message!");
    logger->critical("(4) This is a critical error message!");
  }

  SECTION("test log level") {

    // set the default log level
    logging::set_default_level(logging::level::critical);

    // and make sure we can log with this created object
    logging::info("(5) This should NOT be printed!");
    logging::warn("(5) This should NOT be printed!");
    logging::debug("(5) This should NOT be printed!");
    logging::error("(5) This should NOT be printed!");
    logging::critical("(5) This SHOULD BE printed!!");

    // get a reference to an unknown logger
    auto logger = get_logger("loggerD");

    // now set the default log level for this logger
    logger->set_level(logging::level::critical);

    // now try the various logging functions
    logger->info("(6) This should NOT be printed!");
    logger->warn("(6) This should NOT be printed!");
    logger->debug("(6) This should NOT be printed!");
    logger->error("(6) This should NOT be printed!");
    logger->critical("(6) This SHOULD BE printed!!");

    // and reset it for the next tests
    logging::set_default_level(logging::level::debug);
    logger->set_level(logging::level::critical);
  }

  SECTION("test macro style logging") {

    // these print with the "corsika" logger
    CORSIKA_LOG_INFO("(7) test macro style logging");
    CORSIKA_LOG_DEBUG("(7) test macro style logging");
    CORSIKA_LOG_ERROR("(7) test macro style logging");
    CORSIKA_LOG_CRITICAL("(7) test macro style logging");

    // get a reference to an unknown logger
    auto logger = get_logger("loggerE");

    // add the filename, source information to this logger
    logging::add_source_info(logger);

    // these print with the "loggerE" logger
    CORSIKA_LOGGER_INFO(logger, "(8) test macro style logging");
    CORSIKA_LOGGER_WARN(logger, "(8) test macro {} logging", "style");

    // reset the logging pattern
    logging::reset_pattern(logger);

    // these trace macros should not print file, function, and line
    CORSIKA_LOGGER_TRACE(logger, "(9) test macro style logging:");
  }
}
