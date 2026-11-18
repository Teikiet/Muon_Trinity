/*
 * (c) Copyright 2021 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <corsika/output/BaseOutput.hpp>

namespace corsika {

  /**
   * Generic class to switch off any output.
   *
   * The 'write' method is catch-all and does nothing.
   */

  class WriterOff : public BaseOutput {

  public:
    /**
     * The purpose of this catch-all constructor is to discard all parameters.
     * @tparam TArgs
     */
    template <typename... TArgs>
    WriterOff(TArgs&&...) {}

    virtual ~WriterOff() {}

    void startOfLibrary(boost::filesystem::path const&) override {}

    void endOfShower(unsigned int const) override {}

    void endOfLibrary() override {}

    /**
     * The purpose of this catch-all method is to discard all data.
     * @tparam TArgs
     */
    template <typename... TArgs>
    void write(TArgs&&...) {}

    virtual YAML::Node getConfig() const override { return YAML::Node(); }

  }; // class WriterOff

} // namespace corsika
