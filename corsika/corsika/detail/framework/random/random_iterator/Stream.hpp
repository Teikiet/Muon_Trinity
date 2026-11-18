/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/*
 * Stream.hpp
 *
 *  Created on: 23/02/2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

#include <stdint.h>

#include "detail/tbb/iterators.h"
#include "detail/Engine.hpp"
#include "detail/functors/EngineCaller.hpp"

namespace random_iterator {

  /// philox engine
  typedef detail::philox philox;
  /// threefry  engine
  typedef detail::threefry threefry;
  /// squares3_128  engine
  typedef detail::squares3_128 squares3_128;
  /// squares4_128  engine
  typedef detail::squares4_128 squares4_128;
  /// squares3_64  engine
  typedef detail::squares3_64 squares3_64;
  /// squares4_64  engine
  typedef detail::squares4_64 squares4_64;

#if RANDOM_ITERATOR_R123_USE_AES_NI
  /// ars  engine
  typedef detail::ars ars;
#endif

  /**
   * \class Stream
   *
   * \brief representation for streams of pseudorandom numbers
   *
   * Objects of this class represent collections of up to \f$ {2}^{32} \f$  streams of
   * random numbers for each seed. Each of such streams has a length of \f$ {2}^{64} \f$.
   * If the chosen distribution generates data with 64bits, like ``double``, then each
   * stream can handle up to 128 ExaByte of data. Stream objects are thread-safe and
   * lazy-evaluated. Output is produced only when requested. It is user responsibility to
   * retrieve and store, if needed, the generated output.
   *
   * \tparam DistibutionType STL compliant distribution type. For example:
   * ```std::uniform_real_distribution```, ```std::exponential_distribution```.
   *
   * \tparam EngineType RNG type.
   *
   * There are seven options:
   *
   *  * random_iterator::philox
   *  * random_iterator::threefry
   *  * random_iterator::ars
   *  * random_iterator::squares3_128
   *  * random_iterator::squares4_128
   *  * random_iterator::squares3_64
   *  * random_iterator::squares4_64
   *
   *  Usage:
   *
   * 	1) Direct class instantiation:
   *
   *  \code
   *  \\this code will print the 10 doubles randomly distributed in the range [0.0, 1.0].
   *
   *  \\instantiate real uniform distribution from C++ standard libray
   *  std::uniform_real_distribution<double> udist(0.0, 1.0);
   *
   *  \\
   *  Stream<std::uniform_real_distribution<double>, random_iterator::philox>
   * uniform_stream(udist, 0x1a2b3c4d, 1 );
   *
   *  for( size_t i=0; i<10; ++i )
   *   std::cout << uniform_stream[i] << std::endl;
   *  \endcode
   *
   * 	2) Using ```random_iterator::make_stream(...)```:
   *
   * 	\code
   *  \\this code will print the 10 doubles randomly distributed in the range [0.0, 1.0].
   *
   *  \\instantiate real uniform distribution from C++ standard libray
   *  std::uniform_real_distribution<double> udist(0.0, 1.0);
   *
   *  \\instantiate the generator
   *  random_iterator::philox   rng(0x1a2b3c4d);
   *
   *  \\create stream
   *  auto uniform_stream = random_iterator::make_stream(udist, rng, 1 );
   *
   *  for( size_t i=0; i<10; ++i )
   *   std::cout << uniform_stream[i] << std::endl;
   * 	\endcode
   */
  template <typename DistibutionType, typename EngineType>
  class Stream {
    typedef detail::EngineCaller<DistibutionType, EngineType> caller_type; //!
    typedef random_iterator_tbb::counting_iterator<typename EngineType::advance_type>
        counting_iterator; //!

  public:
    typedef EngineType engine_type;                          //! Engine type
    typedef typename engine_type::state_type state_type;     //! Type of RNG state
    typedef typename engine_type::seed_type seed_type;       //! Type of RNG seed
    typedef typename engine_type::advance_type advance_type; //! Type of RNG displacement
    typedef typename engine_type::init_type init_type;       //! Type of RNG initializer

    typedef DistibutionType distribution_type; //! type of the STL compliant distribution
    typedef typename distribution_type::result_type result_type; //! type of the result

    typedef random_iterator_tbb::transform_iterator<caller_type, counting_iterator>
        iterator; //! type of stream iterator

    Stream() = delete;

    /**
     * Constructor.
     *
     * @param distribution Distribution with stl compliant interface.
     * @param seed Seed for pseudorandom number generation.
     * @param stream Pseudorandom number stream.
     */
    Stream(distribution_type const& distribution, seed_type seed, uint32_t stream = 0)
        : distribution_(distribution)
        , engine_(seed, stream)
        , seed_(seed)
        , stream_(stream) {}

    /**
     * Copy constructor
     *
     * @param other
     */
    Stream(Stream<DistibutionType, EngineType> const& other)
        : distribution_(other.getDistribution())
        , engine_(other.getEngine())
        , seed_(other.getSeed())
        , stream_(other.getStream()) {}

    /**
     * Returns an iterator to the first element of the stream.
     * @return  iterator to the first element of the stream.
     */
    inline iterator begin() const {

      counting_iterator first(0);
      caller_type caller(distribution_, seed_, stream_);

      return iterator(first, caller);
    }

    /**
     * Returns an iterator to the past last element of the stream.
     * @return iterator to the past last element of the stream.
     */
    inline iterator end() const {

      counting_iterator last(std::numeric_limits<advance_type>::max());
      caller_type caller(distribution_, seed_, stream_);

      return iterator(last, caller);
    }

    /**
     * Returns the element at specified location n. No bounds checking is performed.
     *
     * @param n position in the sequence.
     * @return element at location n
     */
    inline result_type operator[](advance_type n) const {

      caller_type caller(distribution_, seed_, stream_);

      return caller(n);
    }

    /**
     *  At each call, this function will produce a pseudorandom number and advance the
     * engine state.
     *
     * @return pseudorandom number distributed according with `DistributionType`
     */
    inline result_type operator()(void) { return distribution_(engine_); }

    /**
     * Get the distribution
     *
     * @return Copy of the object's distribution.
     */
    inline distribution_type getDistribution() const { return distribution_; }

    /**
     * Set the distribution
     * @param dist
     *
     */
    inline void setDistribution(distribution_type dist) { distribution_ = dist; }

    /**
     * Get the associated engine
     *
     * @return Copy of the object's engine
     */
    inline engine_type getEngine() const { return engine_; }

    /**
     * Set the object's engine
     *
     * @param engine
     */
    inline void setEngine(engine_type engine) { engine_ = engine; }

    /**
     * Get the seed
     *
     * @return Copy of the object's seed
     */
    inline const seed_type& getSeed() const { return seed_; }

    /**
     * Set the seed
     *
     * @param seed
     */
    inline void setSeed(const seed_type& seed) { seed_ = seed; }

    /**
     * Get the object's stream number
     *
     * @return stream number
     */
    inline uint32_t getStream() const { return stream_; }

    /**
     * Set the object's stream number
     *
     * @param stream
     */
    inline void setStream(uint32_t stream) { stream_ = stream; }

    friend inline std::ostream& operator<<(
        std::ostream& os, const Stream<DistibutionType, EngineType>& be) {
      return os << "engine: " << be.getEngine() << " stream: " << be.getStream()
                << " distribution: " << be.getDistribution();
    }

  private:
    distribution_type distribution_;
    engine_type engine_;
    seed_type seed_;
    uint32_t stream_;
  };

  /**
   * Stream
   *
   * @param dist STL compliant distribution instance.
   * @param eng Engine instance
   * @param stream Stream number, in the range [0, 2^{32} - 1]
   * @return Stream object.
   */
  template <typename DistibutionType, typename EngineType>
  Stream<DistibutionType, EngineType> make_stream(DistibutionType const& dist,
                                                  EngineType eng, uint32_t stream) {

    return Stream<DistibutionType, EngineType>(dist, eng.getSeed(), stream);
  }

} // namespace random_iterator
