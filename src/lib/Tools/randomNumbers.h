/*!
        @file    randomNumbers.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef RANDOMNUMBERS_INCLUDED
#define RANDOMNUMBERS_INCLUDED

#include "Field/field.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! Base class of random number generators.

/*!
   This class defines the interface of random number
   generator, and implements common methods.
   Practical methods to generate random numbers are
   defined in subclasses.
   This class also implements Gaussian random number and
   method to set a global field of Gaussian random numbers
   and cut it out to the local field for the own node
   (gauss_lex_global()) which is useful in HMC etc.
                                 [25 Dec 2011 H.Matsufuru]
   U1,Z2 are added               [11 Jan 2017 Y.Namekawa]
   Factory is introduced         [ 2 Feb 2017 Y.Namekawa]
 */

class RandomNumbers
{
 public:
  static const std::string class_name;

 protected:
  Bridge::VerboseLevel m_vl;

 public:

  RandomNumbers()
    : m_vl(CommonParameters::Vlevel()) {}

  virtual ~RandomNumbers() {}

 private:
  // non-copyable
  RandomNumbers(const RandomNumbers&);
  RandomNumbers& operator=(const RandomNumbers&);

 public:
  void set_parameter_verboselevel(const Bridge::VerboseLevel vl) { m_vl = vl; }

  virtual double get() = 0;

  void gauss(double& rand1, double& rand2);

  virtual void lex_global(const std::string&, Field&);

  //! gaussian random number defined on global lattice.
  virtual void gauss_lex_global(Field&);

  //! uniform random number defined on global lattice.
  virtual void uniform_lex_global(Field&);

  //! U(1) random number defined on global lattice.
  virtual void U1_lex_global(Field&);

  //! Z(2) random number defined on global lattice.
  virtual void Z2_lex_global(Field&);


  //! gaussian noise for even-odd perconditioned field (S.UEDA)
  virtual void gauss_eo_global(Field&);


  //! save and load random number status.
  virtual void read_file(const std::string&)  = 0;
  virtual void write_file(const std::string&) = 0;


  //! write and read states to/from ofstream (I.Kanamori, for parallel version)
  virtual void read_state(std::ifstream&) {
    vout.crucial(m_vl, "%s: read_state() is not implemented\n", class_name.c_str());
    abort();
  }
  virtual void write_state(std::ofstream&) {
    vout.crucial(m_vl, "%s: write_state() is not implemented\n", class_name.c_str());
    abort();
  }
  virtual void read_state_binary(std::ifstream&) {
    vout.crucial(m_vl, "%s: read_state_binary() is not implemented\n", class_name.c_str());
    abort();
  }
  virtual void write_state_binary(std::ofstream&) {
    vout.crucial(m_vl, "%s: write_state_binary() is not implemented\n", class_name.c_str());
    abort();
  }

  virtual void bcast_state() {
    vout.crucial(m_vl, "%s: bcast_state() is not implemented\n", class_name.c_str());
    abort();
  }

  virtual void pack_state(void *) {
    vout.crucial(m_vl, "%s: pack_state() is not implemented\n", class_name.c_str());
    abort();
  }
  virtual void unpack_state(const void *) {
    vout.crucial(m_vl, "%s: unpack_state() is not implemented\n", class_name.c_str());
    abort();
  }

  //! reset state with new seed.
  virtual void reset(unsigned long seed) = 0;

  //! reset state with new seed and lattice size (I.Kanamaori, for parallel)
  //  virtual void reset(unsigned long seed, std::vector<int>&) = 0;

 protected:
  class rand_gauss_even
  {
   public:
    rand_gauss_even(Field& f, RandomNumbers *rand_gauss)
      : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
    void operator()(const bool do_fill);
    size_t block_size() const;

   private:
    RandomNumbers *m_rand_gauss;
    double *m_ptr;
    size_t m_block;
  };

  class rand_gauss_odd
  {
   public:
    rand_gauss_odd(Field& f, RandomNumbers *rand_gauss)
      : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
    void operator()(const bool do_fill);
    size_t block_size() const;

   private:
    RandomNumbers *m_rand_gauss;
    double *m_ptr;
    size_t m_block;
  };

  //! uniform random number defined on global lattice.
  //void uniform_lex_global_parallel(Field&);

  class rand_uniform
  {
   public:
    rand_uniform(Field& f, RandomNumbers *rand_gauss)
      : m_rand_gauss(rand_gauss), m_ptr(f.ptr(0)), m_block(f.nin()) {}
    void operator()(const bool do_fill);
    size_t block_size() const;

   private:
    RandomNumbers *m_rand_gauss;
    double *m_ptr;
    size_t m_block;
  };

 private:
  template<typename InnerGenerator>
  void generate_global(Field& f);

#ifdef USE_FACTORY
 public:
  typedef RandomNumbers *(*ProductCreator_int)(const int& iseed);
  typedef RandomNumbers *(*ProductCreator_file)(const std::string& filename);

  typedef FactoryTemplate<RandomNumbers, ProductCreator_int>    Factory_int;
  typedef FactoryTemplate<RandomNumbers, ProductCreator_file>   Factory_file;

  static RandomNumbers *New(const IdentifierType& subtype, const int& iseed)
  {
    ProductCreator_int p = Factory_int::Find(subtype);

    return p ? (*p)(iseed) : 0;
  }

  static RandomNumbers *New(const IdentifierType& subtype, const std::string& filename)
  {
    ProductCreator_file p = Factory_file::Find(subtype);

    return p ? (*p)(filename) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
