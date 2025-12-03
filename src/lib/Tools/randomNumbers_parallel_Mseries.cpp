/*!
        @file    randomNumbers_parallel_Mseries.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "randomNumbers.h"
#include "randomNumbers_Mseries.h"
#include "randomNumbers_parallel.h"
#include "randomNumbers_parallel-tmpl.h"


typedef RandomNumbers_Mseries RAND;

template<>
const std::string RandomNumbers_parallel<RandomNumbers_Mseries>::class_name
  = "RandomNumbers_parallel<RancomNumbers_Mseries>";


// use 2-dim parallelizatoin of the generators
template<>
void RandomNumbers_parallel<RAND>::init(const unsigned long s){
  vout.general(m_vl, "%s: init with %lu\n", class_name.c_str(), s);
  init_2d(s);
}

template<>
void RandomNumbers_parallel<RAND>::set_lex_global(Field &v, get_rand &r){
  set_lex_global_2d(v, r);
}

template<>
void RandomNumbers_parallel<RAND>::read_file(const std::string &filename){
  read_file_2d(filename);
}

template<>
void RandomNumbers_parallel<RAND>::write_file(const std::string &filename){
  write_file_2d(filename);
}


// the constructor of Mseries takes int
template<>
RAND* RandomNumbers_parallel<RAND>::new_rand(unsigned long s) const {
  return new RAND((int) s);
}


#ifdef USE_FACTORY
template<>
RandomNumbers * RandomNumbers_parallel<RAND>::create_object_with_int(const int& iseed)
  {
    return new RandomNumbers_parallel<RAND>(iseed);
  }

template<>
RandomNumbers * RandomNumbers_parallel<RAND>::create_object_with_file(const std::string& filename)
  {
    return new RandomNumbers_parallel<RAND>(filename);
  }

template<>
bool RandomNumbers_parallel<RAND>::register_factory()
  {
    bool init1 = RandomNumbers::Factory_int::Register("Mseries_parallel", create_object_with_int);
    bool init2 = RandomNumbers::Factory_file::Register("Mseries_parallel", create_object_with_file);

    return init1 && init2;
  }
#endif


// explicit instanciation
template class
RandomNumbers_parallel<RAND>;

