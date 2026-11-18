/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

#include <map>
#include <cstdint>
#include <random>
#include <string>

#include <corsika/framework/utility/Singleton.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/detail/framework/random/random_iterator/Stream.hpp>
/*!
 * With this class modules can register streams of random numbers.
 */

namespace corsika {

  template <typename CBPRNG = random_iterator::philox>
  class RNGManager : public corsika::Singleton<RNGManager<CBPRNG>> {

    friend class corsika::Singleton<RNGManager<CBPRNG>>;

  public:
    typedef CBPRNG prng_type;
    typedef std::uint64_t seed_type;
    typedef std::string string_type;
    typedef std::map<std::string, prng_type> streams_type;

    RNGManager(RNGManager<CBPRNG> const&) = delete; // since it is a singleton

    RNGManager<CBPRNG>& operator=(RNGManager<CBPRNG> const&) = delete;

    /*!
     * This function is to be called by a module requiring a random-number
     * stream during its initialization.
     *
     * \throws sth. when stream \a pModuleName is already registered
     */
    inline void registerRandomStream(string_type const& streamName);

    /*!
     * returns the pre-stored stream of given name \a pStreamName if
     * available
     */
    inline prng_type& getRandomStream(string_type const& streamName);

    /*!
     * Check whether a stream has been registered.
     */
    inline bool isRegistered(string_type const& streamName) const;

    /*!
     * dumps the names and states of all registered random-number streams
     * into a std::stringstream.
     */
    inline std::stringstream dumpState() const;

    /**
     * @fn const streams_type getRngs&()const
     * @brief Constant access to the streams.
     *
     * @pre
     * @post
     * @return RNGManager::streams_type
     */
    const streams_type& getRngs() const { return rngs_; }

    /**
     * @fn streams_type Rngs()
     * @brief Non-constant access to the streams.
     *
     * @pre
     * @post
     * @return RNGManager::streams_type&
     */
    streams_type& Rngs() { return rngs_; }

    seed_type getSeed() const { return seed_; }

    void setSeed(seed_type seed) {
      seed_ = seed;
      // update the rgn states
      for (auto& [streamName, rng] : rngs_) rng.setSeed(seed_);
    }

  protected:
    RNGManager() = default;

  private:
    streams_type rngs_;
    seed_type seed_;
  };

  typedef typename RNGManager<>::prng_type default_prng_type;

} // namespace corsika

#include <corsika/detail/framework/random/RNGManager.inl>
