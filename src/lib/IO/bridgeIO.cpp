/*!
        @file    bridgeIO.cpp

        @brief

        @author  Satoru Ueda (maintained by I.Kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-05-27 18:31:04 #$

        @version $LastChangedRevision: 2644 $
 */

#include "bridgeIO.h"
#include "Parameters/commonParameters.h"
#include "ResourceManager/threadManager.h"

int Bridge::BridgeIO::m_indent_level = 0;

//====================================================================
// verbose output for c style
// default verbose level, node 0

namespace Bridge {
  //====================================================================
  const std::string BridgeIO::class_name = "BridgeIO";

  //====================================================================
  BridgeIO::BridgeIO(const std::string& filename)
  {
    os_ = NULL;
    m_node_write = -1;
    init(filename);
  }

  //====================================================================
  BridgeIO::BridgeIO(const std::string& filename, int node)
  {
    os_ = NULL;
    m_node_write = -1;
    init(filename, node);
  }

    //====================================================================
  BridgeIO::BridgeIO(const std::string& filename, int node,
                     std::ios_base::openmode mode)
  {
    os_ = NULL;
    m_node_write = -1;
    init(filename, node, mode);
  }


  //====================================================================
  BridgeIO::~BridgeIO()
  {
    tidyup_();
  }


  //====================================================================
  void BridgeIO::init(const std::string& filename)
  {
     // open in all ranks, the default output rank is 0
    int this_node  = Communicator::self();
    init(filename, this_node, std::ios::out);
    m_node_write = 0;  // force the default output rank to 0
  }


  //====================================================================
  void BridgeIO::init(const std::string& filename, int node)
  {
    // open in specified rank, the default output rank is the same
    init(filename, node, std::ios::out);
  }

  //====================================================================
  void BridgeIO::init(const std::string& filename, int node,
                      std::ios_base::openmode mode)
  {
    ThreadManager::assert_single_thread(class_name);
    int this_node  = Communicator::self();

    push_();

    m_node_write = node;

    if( this_node != node) {
      os_ = NULL;
      return;
    }

    if (filename == "stdout") {
      os_ = new std::ostream(std::cout.rdbuf());
    } else {
      os_ = new std::ofstream(filename.c_str(), mode);
    }

    if (!os_) {
      fprintf(stderr, "%s: init: unable to open log file \"%s\".\n", class_name.c_str(), filename.c_str());
      exit(EXIT_FAILURE);

      rewind_();
    }

  }


  //====================================================================
  void BridgeIO::init(const std::ostream& ost)
  {
    int this_node  = Communicator::self();
    init(ost, this_node);
  }

  //====================================================================
  void BridgeIO::init(const std::ostream& ost, int node)
  {
    int this_node  = Communicator::self();

    push_();
    m_node_write = node;

    if( this_node != m_node_write ) {
      os_ = NULL;
      return;
    }

    os_ = new std::ostream(ost.rdbuf());

    if (!os_) {
      fprintf(stderr, "%s: init: unable to open stream.\n", class_name.c_str());
      exit(EXIT_FAILURE);

      rewind_();
    }
  }


  //====================================================================
  void BridgeIO::unset()
  {
    if (os_) delete os_;

    rewind_();
  }


  //====================================================================
  void BridgeIO::push_()
  {
    if ( os_ ) {
      *os_ << std::flush;
    }
    os_info info;
    info.os = os_;
    info.node_write = m_node_write;
    stack_.push(info);
  }


  //====================================================================
  void BridgeIO::rewind_()
  {
    if (stack_.size() > 0) {
      os_info info = stack_.top();
      os_ = info.os;
      m_node_write = info.node_write;

      stack_.pop();
    } else {
      os_ = NULL;
    }
  }


  //====================================================================
  void BridgeIO::tidyup_()
  {
    if (os_) delete os_;

    while (stack_.size() > 0)
    {
      std::ostream *otmp = stack_.top().os;
      if (otmp) delete otmp;
      stack_.pop();
    }

  }


  //====================================================================
  VerboseLevel
  BridgeIO::set_verbose_level(const std::string& str)
  {
    ThreadManager::assert_single_thread(class_name);

    if ((str == "Crucial") || (str == "crucial") || (str == "CRUCIAL")) return Bridge::CRUCIAL;

    if ((str == "General") || (str == "general") || (str == "GENERAL")) return Bridge::GENERAL;

    if ((str == "Detailed") || (str == "detailed") || (str == "DETAILED")) return Bridge::DETAILED;

    if ((str == "Paranoiac") || (str == "paranoiac") || (str == "PARANOIAC")) return Bridge::PARANOIAC;

    if ((str == "NULL") || (str == "null")) return CommonParameters::Vlevel();

    // safe default
    return Bridge::GENERAL;
  }


  //====================================================================
  std::string
  BridgeIO::get_verbose_level(const VerboseLevel vl)
  {
    ThreadManager::assert_single_thread(class_name);

    switch (vl)
    {
    case Bridge::CRUCIAL:
      return "Crucial";

    case Bridge::GENERAL:
      return "General";

    case Bridge::DETAILED:
      return "Detailed";

    case Bridge::PARANOIAC:
      return "Paranoiac";

    default:
      return "NULL";
    }
  }


  //====================================================================
  void
  BridgeIO::crucial(const char *format, ...)
  {
    VerboseLevel vl = CommonParameters::Vlevel();

    if (vl < Bridge::CRUCIAL) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::CRUCIAL, m_node_write, format, arg);
      va_end(arg);
      if( os_ ) *os_ << std::flush;
    }
  }


  //====================================================================
  void
  BridgeIO::general(const char *format, ...)
  {
    VerboseLevel vl = CommonParameters::Vlevel();

    if (vl < Bridge::GENERAL) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::GENERAL, m_node_write, format, arg);
      va_end(arg);
    }
  }


  //====================================================================
  void
  BridgeIO::detailed(const char *format, ...)
  {
    VerboseLevel vl = CommonParameters::Vlevel();

    if (vl < Bridge::DETAILED) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::DETAILED, m_node_write, format, arg);
      va_end(arg);
    }
  }


  //====================================================================
  void
  BridgeIO::paranoiac(const char *format, ...)
  {
    VerboseLevel vl = CommonParameters::Vlevel();

    if (vl < Bridge::PARANOIAC) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::PARANOIAC, m_node_write, format, arg);
      va_end(arg);
    }
  }


  //====================================================================
  // input verbose level, node 0
  void
  BridgeIO::crucial(VerboseLevel vl, const char *format, ...)
  {
    if (vl < Bridge::CRUCIAL) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::CRUCIAL, m_node_write, format, arg);
      va_end(arg);
      if( m_node_write == Communicator::self() && os_ ) *os_ << std::flush;
    }
  }


  //====================================================================
  void
  BridgeIO::general(VerboseLevel vl, const char *format, ...)
  {
    if (vl < Bridge::GENERAL) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::GENERAL, m_node_write, format, arg);
      va_end(arg);
    }
  }


  //====================================================================
  void
  BridgeIO::detailed(VerboseLevel vl, const char *format, ...)
  {
    if (vl < Bridge::DETAILED) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::DETAILED, m_node_write, format, arg);
      va_end(arg);
    }
  }


  //====================================================================
  void
  BridgeIO::paranoiac(VerboseLevel vl, const char *format, ...)
  {
    if (vl < Bridge::PARANOIAC) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      va_start(arg, format);
      print(vl, Bridge::PARANOIAC, m_node_write, format, arg);
      va_end(arg);
    }
  }

  //====================================================================
  // input verbose level, input node
  void
  BridgeIO::crucial(VerboseLevel vl, int node, const char *format, ...)
  {
    if (vl < Bridge::CRUCIAL) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0 ) {
      if (node == m_node_write ) {
        va_start(arg, format);
        print(vl, Bridge::CRUCIAL, node, format, arg);
        va_end(arg);
        if( (m_node_write == Communicator::self()) && os_ ) *os_ << std::flush;
      } else {
#ifdef DEBUG
        if (!isOpen() ){
          // if the output is stdout, it is always open
          // and m_node_rank != node can be by intention
          fprintf(stderr,
                  "BridgeIO::curical: the stream is open for rank %d but specified %d\n",
                  m_node_write, node);
        }
#endif
      }
    }
  }


  //====================================================================
  void
  BridgeIO::general(VerboseLevel vl, int node, const char *format, ...)
  {
    if (vl < Bridge::GENERAL) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      if (node == m_node_write ) {
        va_start(arg, format);
        print(vl, Bridge::GENERAL, node, format, arg);
        va_end(arg);
      } else {
#ifdef DEBUG
        if (!isOpen() ){
          fprintf(stderr,
                  "BridgeIO::general: the stream is open for rank %d but specified %d\n",
                  m_node_write, node);
        }
#endif
      }
    }
  }

  //====================================================================
  void
  BridgeIO::detailed(VerboseLevel vl, int node, const char *format, ...)
  {
    if (vl < Bridge::DETAILED) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      if ( node == m_node_write ) {
        va_start(arg, format);
        print(vl, Bridge::DETAILED, node, format, arg);
        va_end(arg);
      } else {
#ifdef DEBUG
        if (!isOpen() ){
          fprintf(stderr,
                  "BridgeIO::detailed: the stream is open for rank %d but specified %d\n",
                  m_node_write, node);
        }
#endif
      }
    }
  }

  //====================================================================
  void
  BridgeIO::paranoiac(VerboseLevel vl, int node, const char *format, ...)
  {
    if (vl < Bridge::PARANOIAC) return;

    va_list arg;

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      if (node == m_node_write ) {
        va_start(arg, format);
        print(vl, Bridge::PARANOIAC, node, format, arg);
        va_end(arg);
      } else {
#ifdef DEBUG
        if (!isOpen() ){
          fprintf(stderr,
                  "BridgeIO::paranoiac: the stream is open for rank %d but specified %d\n",
                  m_node_write, node);
        }
#endif
      }
    }
  }

  //====================================================================
  std::ostream& BridgeIO::getStream()
  {
    return *os_;
  }


  //====================================================================
  bool BridgeIO::isOpen()
  {
    return os_ && *os_;
  }

  void BridgeIO::increase_indent(){
    int ith = ThreadManager::get_thread_id();
    // synchronization with omp barrier is omitted
    // as m_indent_level is used only in the master thread
    if (ith == 0) {
      ++m_indent_level;
    }
  }

  void BridgeIO::decrease_indent(){
    int ith = ThreadManager::get_thread_id();
    // synchronization with omp barrier is omitted
    // as m_indent_level is used only in the master thread
    if (ith == 0) {
      --m_indent_level;
    }
  }

  void BridgeIO::set_indent(const int level){
    int ith = ThreadManager::get_thread_id();
    // synchronization with omp barrier is omitted
    // as m_indent_level is used only in the master thread
    if (ith == 0) {
      m_indent_level = level;
    }
  }

  //====================================================================
  inline
  void
  BridgeIO::print(VerboseLevel level, VerboseLevel write_level,
                  int node, const char *format, va_list& arg)
  {
    if ((write_level <= level) && (Communicator::nodeid() == node)) {
      if (!os_) {
        std::cerr << "ERROR: BridgeIO: no output stream." << std::endl;
        exit(EXIT_FAILURE);
      }

      vsprintf(buff_, format, arg);

      for (int i = 0; i < m_indent_level; ++i) {
        *os_ << "  ";
      }

      *os_ << buff_;
#ifdef DEBUG
      *os_ << std::flush;
#endif

      if (!os_->good()) {
        std::cerr << "ERROR: BridgeIO: output failed." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }




  //====================================================================
  BridgeIO vout;
}

//====================================================================
//====================================================================
