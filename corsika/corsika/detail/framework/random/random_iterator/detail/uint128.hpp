/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

/*
 * uint128_t.hpp
 *
 *  Created on: 13/03/2021
 *      Author: Antonio Augusto Alves Junior
 */

#pragma once

#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "Macros.h"

namespace random_iterator {

  class uint128_t;
}

namespace std { // This is probably not a good idea
  template <>
  struct is_arithmetic<random_iterator::uint128_t> : std::true_type {};
  template <>
  struct is_integral<random_iterator::uint128_t> : std::true_type {};
  template <>
  struct is_unsigned<random_iterator::uint128_t> : std::true_type {};
} // namespace std

namespace random_iterator {

  // Give uint128_t type traits

  class uint128_t {

  public:
    // Constructors
    uint128_t() = default;

    uint128_t(uint128_t const& rhs) = default;

    uint128_t(uint128_t&& rhs) = default;

    uint128_t(std::string& s);

    uint128_t(const char* s);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    uint128_t(const T& rhs)
#ifdef __BIG_ENDIAN__
        : UPPER(0)
        , LOWER(rhs)
#endif
#ifdef __LITTLE_ENDIAN__
        : LOWER(rhs)
        , UPPER(0)
#endif
    {
      if (std::is_signed<T>::value) {
        if (rhs < T(0)) { UPPER = -1; }
      }
    }

    template <typename S, typename T,
              typename = typename std::enable_if<
                  std::is_integral<S>::value && std::is_integral<T>::value, void>::type>
    uint128_t(const S& upper_rhs, const T& lower_rhs)
#ifdef __BIG_ENDIAN__
        : UPPER(upper_rhs)
        , LOWER(lower_rhs)
#endif
#ifdef __LITTLE_ENDIAN__
        : LOWER(lower_rhs)
        , UPPER(upper_rhs)
#endif
    {
    }

    //  RHS input args only

    // Assignment Operator
    inline uint128_t& operator=(const uint128_t& rhs) = default;

    inline uint128_t& operator=(uint128_t&& rhs) = default;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator=(const T& rhs) {
      UPPER = 0;

      if (std::is_signed<T>::value) {
        if (rhs < T(0)) { UPPER = -1; }
      }

      LOWER = rhs;
      return *this;
    }

    // Typecast Operators
    inline operator bool() const;
    inline operator uint8_t() const;
    inline operator uint16_t() const;
    inline operator uint32_t() const;
    inline operator uint64_t() const;

    // Bitwise Operators
    inline uint128_t operator&(const uint128_t& rhs) const;

    // inline  void export_bits(std::vector<uint8_t> & ret) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator&(const T& rhs) const {
      return uint128_t(0, LOWER & (uint64_t)rhs);
    }

    inline uint128_t& operator&=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator&=(const T& rhs) {
      UPPER = 0;
      LOWER &= rhs;
      return *this;
    }

    inline uint128_t operator|(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator|(const T& rhs) const {
      return uint128_t(UPPER, LOWER | (uint64_t)rhs);
    }

    inline uint128_t& operator|=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator|=(const T& rhs) {
      LOWER |= (uint64_t)rhs;
      return *this;
    }

    inline uint128_t operator^(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator^(const T& rhs) const {
      return uint128_t(UPPER, LOWER ^ (uint64_t)rhs);
    }

    inline uint128_t& operator^=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator^=(const T& rhs) {
      LOWER ^= (uint64_t)rhs;
      return *this;
    }

    inline uint128_t operator~() const;

    // Bit Shift Operators
    inline uint128_t operator<<(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator<<(const T& rhs) const {
      return *this << uint128_t(rhs);
    }

    inline uint128_t& operator<<=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator<<=(const T& rhs) {
      *this = *this << uint128_t(rhs);
      return *this;
    }

    inline uint128_t operator>>(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator>>(const T& rhs) const {
      return *this >> uint128_t(rhs);
    }

    inline uint128_t& operator>>=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator>>=(const T& rhs) {
      *this = *this >> uint128_t(rhs);
      return *this;
    }

    // Logical Operators
    inline bool operator!() const;
    inline bool operator&&(const uint128_t& rhs) const;
    inline bool operator||(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator&&(const T& rhs) const {
      return static_cast<bool>(*this && rhs);
    }

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator||(const T& rhs) const {
      return static_cast<bool>(*this || rhs);
    }

    // Comparison Operators
    inline bool operator==(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator==(const T& rhs) const {
      return (!UPPER && (LOWER == (uint64_t)rhs));
    }

    inline bool operator!=(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator!=(const T& rhs) const {
      return (UPPER | (LOWER != (uint64_t)rhs));
    }

    inline bool operator>(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator>(const T& rhs) const {
      return (UPPER || (LOWER > (uint64_t)rhs));
    }

    inline bool operator<(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator<(const T& rhs) const {
      return (!UPPER) ? (LOWER < (uint64_t)rhs) : false;
    }

    inline bool operator>=(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator>=(const T& rhs) const {
      return ((*this > rhs) | (*this == rhs));
    }

    inline bool operator<=(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline bool operator<=(const T& rhs) const {
      return ((*this < rhs) | (*this == rhs));
    }

    // Arithmetic Operators
    inline uint128_t operator+(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator+(const T& rhs) const {
      return uint128_t(UPPER + ((LOWER + (uint64_t)rhs) < LOWER), LOWER + (uint64_t)rhs);
    }

    inline uint128_t& operator+=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator+=(const T& rhs) {
      return *this += uint128_t(rhs);
    }

    inline uint128_t operator-(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator-(const T& rhs) const {
      return uint128_t((uint64_t)(UPPER - ((LOWER - rhs) > LOWER)),
                       (uint64_t)(LOWER - rhs));
    }

    inline uint128_t& operator-=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator-=(const T& rhs) {
      return *this = *this - uint128_t(rhs);
    }

    inline uint128_t operator*(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator*(const T& rhs) const {
      return *this * uint128_t(rhs);
    }

    inline uint128_t& operator*=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator*=(const T& rhs) {
      return *this = *this * uint128_t(rhs);
    }

    inline uint128_t operator/(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator/(const T& rhs) const {
      return *this / uint128_t(rhs);
    }

    inline uint128_t& operator/=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator/=(const T& rhs) {
      return *this = *this / uint128_t(rhs);
    }

    inline uint128_t operator%(const uint128_t& rhs) const;

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t operator%(const T& rhs) const {
      return *this % uint128_t(rhs);
    }

    inline uint128_t& operator%=(const uint128_t& rhs);

    template <typename T,
              typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
    inline uint128_t& operator%=(const T& rhs) {
      return *this = *this % uint128_t(rhs);
    }

    // Increment Operator
    inline uint128_t& operator++();
    inline uint128_t operator++(int);

    // Decrement Operator
    inline uint128_t& operator--();
    inline uint128_t operator--(int);

    // Nothing done since promotion doesn't work here
    inline uint128_t operator+() const;

    // two's complement
    inline uint128_t operator-() const;

    // Get private values
    inline uint64_t upper() const;
    inline uint64_t lower() const;

    // Get bitsize of value
    inline uint8_t bits() const;

    // rotate right 64
    inline uint128_t rotate_right() const;

    // Get string representation of value
    inline std::string str(uint8_t base = 10, const unsigned int& len = 0) const;

  private:
    std::pair<uint128_t, uint128_t> divmod(const uint128_t& lhs,
                                           const uint128_t& rhs) const;
    void init(const char* s);
    void _init_hex(const char* s);
    void _init_dec(const char* s);
    void _init_oct(const char* s);

#if defined(__powerpc64__) || defined(__x86_64__)
    static inline uint128_t mul128(uint128_t const x, uint128_t const y) {
      uint128_t z;
#ifdef __powerpc64__
      asm("mulhdu %3\n\t" : "=a"(z.LOWER), "=d"(z.UPPER) : "%0"(x.LOWER), "rm"(y.LOWER));
#else
      asm("mulq %3\n\t" : "=a"(z.LOWER), "=d"(z.UPPER) : "%0"(x.LOWER), "rm"(y.LOWER));
#endif
      z.UPPER += (x.UPPER * y.LOWER) + (x.LOWER * y.UPPER);
      return z;
    }
#endif

#ifdef __BIG_ENDIAN__
    uint64_t UPPER, LOWER;
#endif
#ifdef __LITTLE_ENDIAN__
    uint64_t LOWER, UPPER;
#endif
  };

  // lhs type T as first arguemnt
  // If the output is not a bool, casts to type T

  // Bitwise Operators
  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator&(const T& lhs, const uint128_t& rhs) {
    return rhs & lhs;
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator&=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(rhs & lhs);
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator|(const T& lhs, const uint128_t& rhs) {
    return rhs | lhs;
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator|=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(rhs | lhs);
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator^(const T& lhs, const uint128_t& rhs) {
    return rhs ^ lhs;
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator^=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(rhs ^ lhs);
  }

  // Bitshift operators
  inline uint128_t operator<<(const bool& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const uint8_t& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const uint16_t& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const uint32_t& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const uint64_t& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const int8_t& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const int16_t& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const int32_t& lhs, const uint128_t& rhs);
  inline uint128_t operator<<(const int64_t& lhs, const uint128_t& rhs);

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator<<=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(uint128_t(lhs) << rhs);
  }

  inline uint128_t operator>>(const bool& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const uint8_t& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const uint16_t& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const uint32_t& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const uint64_t& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const int8_t& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const int16_t& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const int32_t& lhs, const uint128_t& rhs);
  inline uint128_t operator>>(const int64_t& lhs, const uint128_t& rhs);

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator>>=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(uint128_t(lhs) >> rhs);
  }

  // Comparison Operators
  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline bool operator==(const T& lhs, const uint128_t& rhs) {
    return (!rhs.upper() && ((uint64_t)lhs == rhs.lower()));
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline bool operator!=(const T& lhs, const uint128_t& rhs) {
    return (rhs.upper() | ((uint64_t)lhs != rhs.lower()));
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline bool operator>(const T& lhs, const uint128_t& rhs) {
    return (!rhs.upper()) && ((uint64_t)lhs > rhs.lower());
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline bool operator<(const T& lhs, const uint128_t& rhs) {
    if (rhs.upper()) { return true; }
    return ((uint64_t)lhs < rhs.lower());
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline bool operator>=(const T& lhs, const uint128_t& rhs) {
    if (rhs.upper()) { return false; }
    return ((uint64_t)lhs >= rhs.lower());
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline bool operator<=(const T& lhs, const uint128_t& rhs) {
    if (rhs.upper()) { return true; }
    return ((uint64_t)lhs <= rhs.lower());
  }

  // Arithmetic Operators
  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator+(const T& lhs, const uint128_t& rhs) {
    return rhs + lhs;
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator+=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(rhs + lhs);
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator-(const T& lhs, const uint128_t& rhs) {
    return -(rhs - lhs);
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator-=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(-(rhs - lhs));
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator*(const T& lhs, const uint128_t& rhs) {
    return rhs * lhs;
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator*=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(rhs * lhs);
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator/(const T& lhs, const uint128_t& rhs) {
    return uint128_t(lhs) / rhs;
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator/=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(uint128_t(lhs) / rhs);
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline uint128_t operator%(const T& lhs, const uint128_t& rhs) {
    return uint128_t(lhs) % rhs;
  }

  template <typename T,
            typename = typename std::enable_if<std::is_integral<T>::value, T>::type>
  inline T& operator%=(T& lhs, const uint128_t& rhs) {
    return lhs = static_cast<T>(uint128_t(lhs) % rhs);
  }

  // IO Operator
  inline std::ostream& operator<<(std::ostream& stream, const uint128_t& rhs);

} // namespace random_iterator

#include "uint128.inl"
