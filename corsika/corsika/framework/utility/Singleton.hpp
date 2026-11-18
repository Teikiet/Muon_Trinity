/*
 * (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */

#pragma once

namespace corsika {
  /**
   * \class Singleton Singleton.h utl/Singleton.h
   *
   * \brief Curiously Recurring Template Pattern (CRTP) for Meyers singleton
   *
   * The singleton class is implemented as follows
   * \code
   * #include <utl/Singleton.hpp>
   *
   * class SomeClass : public utl::Singleton<SomeClass> {
   *   ...
   * private:
   *   // prevent creation, destruction
   *   SomeClass() { }
   *   ~SomeClass() { }
   *
   *   friend class utl::Singleton<SomeClass>;
   * };
   * \endcode
   * Singleton automatically prevents copying of the derived class.
   *
   * \author Darko Veberic
   * \date 9 Aug 2006
   * \version $Id: Singleton.h 25091 2014-01-30 09:49:57Z darko $
   * \ingroup stl
   */

  template <typename T>
  class Singleton {
  public:
    static T& getInstance() {
      static T instance;
      return instance;
    }

    Singleton(const Singleton&) =
        delete; // Singleton Classes should not be copied. Removes move constructor and
                // move assignment as well
    Singleton& operator=(const Singleton&) =
        delete; // Singleton Classes should not be copied.

  protected:
    // derived class can call ctor and dtor
    Singleton() {}
    ~Singleton() {}
  };

} // namespace corsika
