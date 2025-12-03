/*!
        @file    randomNumbers_parallel.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef RANDOMNUMBERS_PARALLEL_INCLUDED
#define RANDOMNUMBERS_PARALLEL_INCLUDED

#include <assert.h>
#include <string>

#include "randomNumbers.h"
#include "Communicator/communicator.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Random number generator: parallel version

template<typename RAND>
class RandomNumbers_parallel : public RandomNumbers
{

 public:
  static const std::string class_name;

 private:
  unique_ptr<RAND> m_rand_master;
  std::vector< unique_ptr<RAND> > m_rand_array;
  static constexpr int themalization = 1000; // disturb the internal states
  virtual void init(const unsigned long s);
  virtual void init(const std::vector<unsigned long>& key);

 public:
  RandomNumbers_parallel(const unsigned long s);
  RandomNumbers_parallel(const unsigned long s, const std::vector<int> &local_lattice);
  RandomNumbers_parallel(const std::vector<unsigned long>& key);
  RandomNumbers_parallel(const std::vector<unsigned long>& key, const std::vector<int> &local_lattice);
  RandomNumbers_parallel(const std::string& filename){
    read_file(filename);
  };


  double get() {
    return m_rand_master->get();
  }

  void gauss(double& rand1, double& rand2) {
    m_rand_master->gauss(rand1, rand2);
  }

  double get_block();

  void lex_global(const std::string&, Field&);

  //! gaussian random number defined on global lattice.
  void gauss_lex_global(Field&);

  //! uniform random number defined on global lattice.
  void uniform_lex_global(Field&);

  //! U(1) random number defined on global lattice.
  void U1_lex_global(Field&);

  //! Z(2) random number defined on global lattice.
  void Z2_lex_global(Field&);

  //! gaussian noise for even-odd perconditioned field (S.UEDA)
  void gauss_eo_global(Field&);


  void write_file(const std::string&);
  void read_file(const std::string&);

  void reset(unsigned long seed);
  void reset(unsigned long seed, const std::vector<int> &local_lattice);

protected:

  // interface to create a random number generator
  RAND* new_rand(unsigned long s) const {
    return new RAND(s);
  }
  RAND* new_rand(const std::vector<unsigned long>&) const {
    return nullptr;
  }

  std::vector<int> m_Nsize;

  void init_2d(unsigned long seed);


  // interface to obtain random numbers
  class get_rand {
  public:
    get_rand(): m_block_size(1) { };
    get_rand(int size): m_block_size(size) { };

    virtual double get(RandomNumbers &generator)=0;
    virtual void get_block(double *rands, RandomNumbers &generator){
      for(int i=0; i<m_block_size; ++i){
        rands[i]=get(generator);
      }
    }
  protected:
    int m_block_size;
  };

  class get_gauss : public get_rand {
  public:
    get_gauss(int size) : get_rand(size), need_generate(true), r_stored(0) { }
    double get(RandomNumbers &generator){
      if(!need_generate){
        need_generate=true;
        return r_stored;
      }
      double r1, r2;
      generator.gauss(r1,r2);
      r_stored=r2;
      need_generate=false;
      return r1;
    }
  private:
    bool need_generate;
    double r_stored;
  };

  class get_uniform : public get_rand {
  public:
    double get(RandomNumbers &generator){
      return generator.get();
    }
  };

  class get_U1: public get_rand {
  public:
    get_U1() : get_rand(2) { }

    get_U1(int size) : get_rand(size) {
      if(size %2 !=0){
        vout.crucial("get_u1: the size must be even (given: %d)\n", size);
        exit(EXIT_FAILURE);
      }
    }
    double get(RandomNumbers &generator){
      vout.crucial("get_u1::get() is not supported, try get_u1::get_block() instead.\n");
      exit(EXIT_FAILURE);
    }
    void get_block(double *rands, RandomNumbers &generator){
      for(int i=0; i<get_rand::m_block_size/2; ++i){
        double r1, r2;
        generator.gauss(r1,r2);
        double arg=atan2(r2,r1);
        rands[2*i]   = cos(arg);
        rands[2*i+1] = sin(arg);
      }
    }

  };

  class get_Z2 : public get_rand {
  public:
    get_Z2(int size) : get_rand(size) {}
    double get(RandomNumbers &generator){
      double rn1 = generator.get();
      double rn2 = floor(2.0 * rn1);
      return (2.0 * rn2 - 1.0);
    }
  };

  void set_lex_global(Field&, get_rand&);
  //  void write_file(const std::string& filename);
  //  void read_file(const std::string& filename);

protected:
  void set_lex_global_2d(Field&, get_rand&);
  void write_file_2d(const std::string& filename);
  void read_file_2d(const std::string& filename);

  // to be implemented
  void set_lex_global_4d(Field&, get_rand&);
  void write_file_4d(const std::string& filename);
  void read_file_4d(const std::string& filename);


#ifdef USE_FACTORY
 private:
  static RandomNumbers *create_object_with_int(const int& iseed);
  static RandomNumbers *create_object_with_file(const std::string& filename);

 public:
  static bool register_factory();

#endif

};


#endif
