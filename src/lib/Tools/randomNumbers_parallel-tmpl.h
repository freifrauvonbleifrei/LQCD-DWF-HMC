/*!
        @file    randomNumbers_parallel-tmpl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef RANDOMNUMBERS_PARALLEL_TMPL_INCLUDED
#define RANDOMNUMBERS_PARALLEL_TMPL_INCLUDED

#include "lib/Tools/timer.h"

/*
namespace {
  class get_rand(){
    get_rand(): m_block_size(1);
    virtual get_rand(int size): m_block_size(size);
    virtual double get(RandomNumbers &generator);
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
    get_gauss(int size) : get_rand(size), need_genarate(true), r_stored(0);
    double get(RandomNumbers &generator){
      if(!need_generate){
        need_generate=true;
        return r_stored;
      }

      double r1, r2;
      generator.rand(r1,r2);
      r_stored=r2;
      need_generate=false;
      return r1;
    }
  private:
    bool need_generate;
    double r_strored;
  };

class get_uniform : public get_rand {
public:
  double get(RandomNumbers &generator){
      retrun generator.get();
    }
  };

class get_u1: public get_rand {
  public:
    get_u1(int size) : get_rand(size), need_genarate(true), r_stored(0){
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
      for(int i=0; i<m_block_size/2; ++i){
        double r1=generator.get();
        double r2=generator.get();
        double arg=atan2(r2,r1);
        rands[2*i]   = cos(arg);
        rands[2*i+1] = sin(arg);
      }
    }


  };

class get_z2 : public get_rand {
public:
  double get(RandomNumbers &generator){
    double rn1 = generator.get();
    double rn2 = floor(2.0 * rn1);
    return (2.0 * rn2 - 1.0);
  }
};

}
*/


template<typename RAND>
RandomNumbers_parallel<RAND>::RandomNumbers_parallel(const unsigned long s) {
  m_Nsize.resize(4);
  m_Nsize[0] = CommonParameters::Nx();
  m_Nsize[1] = CommonParameters::Ny();
  m_Nsize[2] = CommonParameters::Nz();
  m_Nsize[3] = CommonParameters::Nt();
  init(s);
};

template<typename RAND>
RandomNumbers_parallel<RAND>::RandomNumbers_parallel(const unsigned long s, const std::vector<int> &local_lattice){
  m_Nsize.resize(4);
  m_Nsize[0] = local_lattice[0];
  m_Nsize[1] = local_lattice[1];
  m_Nsize[2] = local_lattice[2];
  m_Nsize[3] = local_lattice[3];
  init(s);
}

template<typename RAND>
RandomNumbers_parallel<RAND>::RandomNumbers_parallel(const std::vector<unsigned long>& key) {
  m_Nsize[0] = CommonParameters::Nx();
  m_Nsize[1] = CommonParameters::Ny();
  m_Nsize[2] = CommonParameters::Nz();
  m_Nsize[3] = CommonParameters::Nt();
  init(key);
};

template<typename RAND>
RandomNumbers_parallel<RAND>::RandomNumbers_parallel(const std::vector<unsigned long>& key, const std::vector<int> &local_lattice){
  m_Nsize.resize(4);
  m_Nsize[0] = local_lattice[0];
  m_Nsize[1] = local_lattice[1];
  m_Nsize[2] = local_lattice[2];
  m_Nsize[3] = local_lattice[3];
  init(key);
}

template<typename RAND>
void RandomNumbers_parallel<RAND>::init(const std::vector<unsigned long>& key) {
  vout.crucial("%s: initialization with key is not implemented.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}


template<typename RAND>
void RandomNumbers_parallel<RAND>::reset(unsigned long seed){
  init(seed);
}

template<typename RAND>
void RandomNumbers_parallel<RAND>::reset(unsigned long seed, const std::vector<int> &local_lattice){
  m_Nsize[0] = local_lattice[0];
  m_Nsize[1] = local_lattice[1];
  m_Nsize[2] = local_lattice[2];
  m_Nsize[3] = local_lattice[3];
  init(seed);
}

template<typename RAND>
void RandomNumbers_parallel<RAND>::init_2d(const unsigned long seed){
  vout.detailed("%s: init_2d is called\n", class_name.c_str());

  const int Nx = m_Nsize[0];
  const int Ny = m_Nsize[1];
  const int Nz = m_Nsize[2];
  const int Nt = m_Nsize[3];

  const int Lx = Nx*Communicator::npe(0);
  const int Ly = Ny*Communicator::npe(1);
  const int Lz = Nz*Communicator::npe(2);
  const int Lt = Nt*Communicator::npe(3);

  const int ipe_x = Communicator::ipe(0);
  const int ipe_y = Communicator::ipe(1);
  const int ipe_z = Communicator::ipe(2);
  const int ipe_t = Communicator::ipe(3);

  if(m_rand_array.size() != Nz*Nt){
    m_rand_array.resize(Nz*Nt);
  }
  m_rand_master.reset(new_rand(seed));

  // supress the message
  Bridge::VerboseLevel vlevel_keep = CommonParameters::Vlevel();
  CommonParameters::init_Vlevel(Bridge::VerboseLevel::CRUCIAL);
#pragma omp parallel for collapse(2)
  for(int t=0; t<Nt; ++t){
    for(int z=0; z<Nz; ++z){
      int il_z = ipe_z*Nz + z;
      int il_t = ipe_t*Nt + t;
      int il_zt= il_z + Lz*il_t;
      int zt = z + Nz*t;
      unsigned long seed_zt = (seed+1)*Lz*Lt + il_zt;
      m_rand_array[zt].reset(new_rand(seed_zt));
      for(int i=0; i<themalization; ++i){
        m_rand_array[zt]->get();
      }
    }
  }
  CommonParameters::init_Vlevel(vlevel_keep);
}

template<typename RAND>
void RandomNumbers_parallel<RAND>::set_lex_global_2d(Field &v, get_rand& distribution){
    const int Nx = m_Nsize[0];
    const int Ny = m_Nsize[1];
    const int Nz = m_Nsize[2];
    const int Nt = m_Nsize[3];

    const int Lx = Nx*Communicator::npe(0);
    const int Ly = Ny*Communicator::npe(1);
    const int Lz = Nz*Communicator::npe(2);
    const int Lt = Nt*Communicator::npe(3);

    const int ipe_x = Communicator::ipe(0);
    const int ipe_y = Communicator::ipe(1);
    const int ipe_z = Communicator::ipe(2);
    const int ipe_t = Communicator::ipe(3);

    const int Nex=v.nex();
    const int Nin=v.nin();

    std::vector<double> rand_block(Nin);

    for (int j = 0; j < Nex; ++j) {
      bool in_j = true;

#pragma omp for collapse(2)
      for(int t=0; t<Nt; ++t){
        for(int z=0; z<Nz; ++z){
          int zt = z + Nz*t;
          int site = Nx*Ny*(z + Nz*t);
          for(int y=0; y<Ly; ++y){
            bool in_y = in_j && (y >= ipe_y * Ny) && (y < (ipe_y + 1) * Ny);
            for(int x=0; x<Lx; ++x){
              bool in_x = in_y && (x >= ipe_x * Nx) && (x < (ipe_x + 1) * Nx);

              if(in_x){
                assert(site<Nx*Ny*Nz*Nt);
                distribution.get_block(&rand_block[0], *m_rand_array[zt]);
                for(int i=0; i<Nin; ++i){
                  double r = rand_block[i];
                  v.set(i,   site, j, r);
                }
                site++;
              } else {
                distribution.get_block(&rand_block[0], *m_rand_array[zt]);
              }
            }} // x,y
        }}// z,t (omp parallel for)
     }// j
}



//====================================================================
template<typename RAND>
void RandomNumbers_parallel<RAND>::lex_global(const std::string& str_rand_type, Field& f)
{
  if (str_rand_type == "Gaussian") {
    gauss_lex_global(f);
  } else if (str_rand_type == "Uniform") {
    uniform_lex_global(f);
  } else if (str_rand_type == "U1") {
    U1_lex_global(f);
  } else if (str_rand_type == "Z2") {
    Z2_lex_global(f);
  } else {
    vout.crucial(m_vl, "Error at %s: unsupported rand_type \"%s\"\n", class_name.c_str(), str_rand_type.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename RAND>
void RandomNumbers_parallel<RAND>::gauss_lex_global(Field &v){
  int nin=v.nin();
#pragma omp parallel
  {
    get_gauss distribution(nin);
    set_lex_global(v, distribution);
  }
}

//====================================================================
template<typename RAND>
void RandomNumbers_parallel<RAND>::uniform_lex_global(Field &v){
#pragma omp parallel
  {
    get_uniform distribution;
    set_lex_global(v, distribution);
  }
}

//====================================================================
template<typename RAND>
void RandomNumbers_parallel<RAND>::U1_lex_global(Field &v){
  int nin=v.nin();
#pragma omp parallel
  {
    get_U1 distribution(nin);
    set_lex_global(v, distribution);
  }
}

//====================================================================
template<typename RAND>
void RandomNumbers_parallel<RAND>::Z2_lex_global(Field &v){
  int nin=v.nin();

  // assumes the field is complex
  assert(nin % 2 == 0);

#pragma omp parallel
  {
    get_Z2 distribution(nin);
    set_lex_global(v, distribution);
#pragma omp barrier
    double RF2 = 1.0 / sqrt(2.0);
    scal(v, RF2);
  }
}

//====================================================================
template<typename RAND>
void RandomNumbers_parallel<RAND>::gauss_eo_global(Field &v){
  vout.crucial("%s: gauss_eo_global is not implemented.\n", class_name.c_str());
  exit(EXIT_FAILURE);
}

template<typename RAND>
void RandomNumbers_parallel<RAND>::read_file_2d(const std::string &filename){
  vout.general(m_vl, "%s: read from file = %s\n", class_name.c_str(), filename.c_str());

  Timer timer("reading random field");
  timer.start();

  const int Nz = m_Nsize[2];
  const int Nt = m_Nsize[3];

  const int ipe_z = Communicator::ipe(2);
  const int ipe_t = Communicator::ipe(3);

  const int npe_x = Communicator::npe(0);
  const int npe_y = Communicator::npe(1);
  const int npe_z = Communicator::npe(2);
  const int npe_t = Communicator::npe(3);

  const int myrank=Communicator::self();

  // open file
  std::ifstream infile;
  if (Communicator::is_primary()) {
    infile.open(filename.c_str(), std::ios::in | std::ios::binary);

    if (!infile) {
      vout.crucial(m_vl, "Error at %s: unable to open input file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }

  }


  // serial random number
  if (Communicator::is_primary()) {
    m_rand_master->read_state_binary(infile);
  }
  m_rand_master->bcast_state();

  // parallel random numbers
  if(Communicator::size() == 1){
    for(int i=0; i<m_rand_array.size(); ++i){
      m_rand_array[i]->read_state_binary(infile);
    }
  } else {
    RAND rand_tmp(0);
    std::vector< typename RAND::state_type > buffer(Nz);
    size_t buffer_size=sizeof(typename RAND::state_type)*Nz;
    const int recv_from=0; // primary
    int send_to;
    for(int pt=0; pt<npe_t; ++pt){
      for(int t=0; t<Nt; ++t){
        //        int gt = ipe_t*Nt + t;
        for(int pz=0; pz<npe_z; ++pz){
          if (Communicator::is_primary()) {
            for(int z=0; z<Nz; ++z){
              rand_tmp.read_state_binary(infile);
              rand_tmp.pack_state(&buffer[z]);
            }
          }
          for(int px=0; px<npe_x; ++px){
            for(int py=0; py<npe_y; ++py){
              int grid_coord[4] = {px, py, pz, pt};
              Communicator::grid_rank(&send_to, grid_coord);
              int tag = send_to;
              if(send_to != 0){
                Communicator::Base::send_1to1(buffer_size, &buffer[0], &buffer[0], send_to, recv_from, tag);
              }
              if(myrank == send_to){
                for(int z=0; z<Nz; ++z){
                  int zt = z + Nz*t;
                  m_rand_array[zt]->unpack_state(&buffer[z]);
                }
              }
            }
          }
        } // pt
      } // z
    } // pz
  }
  if (Communicator::is_primary()) {
    infile.close();
  }
  Communicator::sync();
  timer.stop();
  timer.report();
}

template<typename RAND>
void RandomNumbers_parallel<RAND>::write_file_2d(const std::string &filename){

  vout.general(m_vl, "%s: write down to file = %s\n", class_name.c_str(), filename.c_str());

  Timer timer("writing random field");
  timer.start();

  const int Nz = m_Nsize[2];
  const int Nt = m_Nsize[3];
  const int ipe_z = Communicator::ipe(2);
  const int ipe_t = Communicator::ipe(3);
  const int npe_z = Communicator::npe(2);
  const int npe_t = Communicator::npe(3);
  const int myrank=Communicator::self();

  // open file
  std::ofstream outfile;
  if (Communicator::is_primary()) {
    outfile.open(filename.c_str(), std::ios::out | std::ios::binary);
    if (!outfile) {
      vout.crucial(m_vl, "Error at %s: unable to open output file.\n", class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }

  // serial random number
  if (Communicator::is_primary()) {
    m_rand_master->write_state_binary(outfile);
  }

  // parallel random numbers
  if(Communicator::size() == 1){
    for(int i=0; i<m_rand_array.size(); ++i){
      m_rand_array[i]->write_state_binary(outfile);
    }
  } else {
    RAND rand_tmp(0);
    std::vector< typename RAND::state_type > buffer(Nz);
    size_t buffer_size=sizeof(typename RAND::state_type)*Nz;
    int recv_from;
    const int send_to=0; // primary
    for(int pt=0; pt<npe_t; ++pt){
      for(int t=0; t<Nt; ++t){
        for(int pz=0; pz<npe_z; ++pz){
          int grid_coord[4] = {0, 0, pz, pt};
          Communicator::grid_rank(&recv_from, grid_coord);
          if(myrank == recv_from){
            for(int z=0; z<Nz; ++z){
              int zt = z + Nz*t;
              m_rand_array[zt]->pack_state(&buffer[z]);
            }
          }
          int tag = recv_from;
          if(recv_from != 0){
            Communicator::Base::send_1to1(buffer_size, &buffer[0], &buffer[0], send_to, recv_from, tag);
          }
          if (Communicator::is_primary()) {
            for(int z=0; z<Nz; ++z){
              rand_tmp.unpack_state(&buffer[z]);
              rand_tmp.write_state_binary(outfile);
            }
          }
          Communicator::sync();
        } // pt
      } // z
    } // pz
  }
  if (Communicator::is_primary()) {
    outfile.close();
  }
  timer.stop();
  timer.report();
}






#endif
