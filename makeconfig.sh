#!/bin/sh

echo ""
echo "################################################################################"
echo "######################## Bridge++ Makefile Setup script ########################"
echo "################################################################################"
#
echo ""

datestr=`date +%Y%m%d%H%M`
if [ -f Makefile ]; then
    echo "Current Makefile is saved to Makefile.$datestr."
    cp -p Makefile Makefile.$datestr
    echo ""
fi

#----------------------------------------------------------------
echo "### Select Platform ###"
echo "1 [PC WorkStation with GNU C++] (default)"
echo "2 [PC WorkStation with Clang C++]"
echo "3 [PC WorkStation with Intel C++]"
echo "4 [PC WorkStation with NVIDIA HPC SDK]"
echo "5 [Fugaku_CLANG]"
echo "6 [SX_Aurora@KEK]"
echo ""

echo -n "Platform? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	TARGET_FLAG="PC_GNU"
	echo "You have selected [PC Workstation with GNU C++]"
	;;
    2)
	TARGET_FLAG="PC_CLANG"
	echo "You have selected [PC Workstation with Clang C++]"
	;;
    3)
	TARGET_FLAG="PC_INTEL"
	echo "You have selected [PC Workstation with Intel C++]"
	;;
    4)
	TARGET_FLAG="PC_NVIDIA"
	echo "You have selected [PC Workstation with NVIDIA HPC SDK]"
	;;
    5)
	TARGET_FLAG="Fugaku_CLANG"
	echo "You have selected [Fugaku_CLANG]"
	;;
    6)
	TARGET_FLAG="SX_AURORA"
	echo "You have selected [SX_AURORA]"
	;;
    *)
	TARGET_FLAG="PC_GNU"
	echo "You have selected default value [PC Workstation With GNU C++]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Select Communicator Type ###"
echo "1 [Single] (default)"
echo "2 [MPI]"
echo ""

echo -n "Communication Type? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	MPI_FLAG="no"
	echo "You have selected [Single]"
	;;
    2)
	MPI_FLAG="yes"
	echo "You have selected [MPI]"
	;;
    *)
	MPI_FLAG="no"
	echo "You have selected default value [Single]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Select Multi-threading Option ###"
echo "1 [OpenMP]"
echo "2 [No] (default)"
echo ""

echo -n "Use Multi-threading? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	MULTITHREAD_FLAG="omp"
	echo "You have selected [OpenMP]"
	;;
    2)
	MULTITHREAD_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	MULTITHREAD_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Optimized Code ###"

if [ "${MULTITHREAD_FLAG}" != "no" ]; then

    echo "## Optimized code is required when multi-threading is enabled."
    echo
    echo "[Yes] enforced."

    OPT_CODE_FLAG="yes"

else
    echo "1 [Yes] (default)"
    echo "2 [No]"
    echo ""

    echo -n "Use Optimized Code? > "
    read VALUE

    echo ""

    case ${VALUE} in
	1)
	    OPT_CODE_FLAG="yes"
	    echo "You have selected [Yes]"
	    ;;
	2)
	    OPT_CODE_FLAG="no"
	    echo "You have selected [No]"
	    ;;
	*)
	    OPT_CODE_FLAG="yes"
	    echo "You have selected default value [Yes]"
	    ;;
    esac
fi

echo ""

#----------------------------------------------------------------
echo "### Select Compiler Optimization Level ###"
echo "1 [Low]"
echo "2 [Standard] (default)"
echo "3 [High]"
echo "4 [Trace] only for SX-Aurora"
echo ""

echo -n "Optimization Level? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	OPT_LEVEL_FLAG="low"
	echo "You have selected [Low]"
	;;
    2)
	OPT_LEVEL_FLAG="std"
	echo "You have selected [Standard]"
	;;
    3)
	OPT_LEVEL_FLAG="high"
	echo "You have selected [High]"
	;;
    4)
	OPT_LEVEL_FLAG="trace"
	echo "You have selected [Trace]"
	;;
    *)
	OPT_LEVEL_FLAG="std"
	echo "You have selected default value [Standard]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Static linkage Option ###"
echo "## [Yes] needs static system libraries"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Static linkage? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	STATIC_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	STATIC_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	STATIC_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Debug Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Debug Option? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	DEBUG_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	DEBUG_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	DEBUG_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Select Complex Type ###"
echo "1 [STL complex] (default)"
echo "2 [C99 complex]"
echo ""

echo -n "Complex Type? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	COMPLEX_FLAG="std"
	echo "You have selected [STL complex]"
	;;
    2)
	COMPLEX_FLAG="c99"
	echo "You have selected [C99 complex]"
	;;
    *)
	COMPLEX_FLAG="std"
	echo "You have selected default value [STL complex]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Select Use C++11 Functions ###"
if [ $TARGET_FLAG = "FX1000" ]; then
    CPP11_FLAG="no"
    echo ""
    echo "CPP11 does not work in $TARGET_FLAG."
    echo "Automatically select no."
else
    echo "1 [Yes] (default)"
    echo "2 [No]"
    echo ""

    echo -n "Use cpp11? > "
    read VALUE

    echo ""

    case ${VALUE} in
	1)
	    CPP11_FLAG="yes"
	    echo "You have selected [Yes]"
	    ;;
	2)
	    CPP11_FLAG="no"
	    echo "You have selected [No]"
	    ;;
	*)
	    CPP11_FLAG="yes"
	    echo "You have selected default value [Yes]"
	    ;;
    esac
fi
echo ""

#----------------------------------------------------------------
echo "### Select Enable Python Interface ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use python? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	PYTHON_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	PYTHON_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	PYTHON_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo

#----------------------------------------------------------------
echo "### Use Alternative Code ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Alternative Code? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	ALT_CODE_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	ALT_CODE_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	ALT_CODE_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo

#----------------------------------------------------------------
echo "### Build directory ###"
echo
echo "which directory to build and store library and include files?"
echo -n "[default: build] > "
read VALUE

if [ -n "${VALUE}" ]; then
    BUILD_DIR=${VALUE}
else
    BUILD_DIR="build"
fi

echo

#----------------------------------------------------------------
echo "### Build Tests directory ###"
echo
echo "which directory to build Tests and store executables?"
echo -n "[default: ${BUILD_DIR}/tests] > "
read VALUE

if [ -n "${VALUE}" ]; then
    TEST_DIR=${VALUE}
else
    TEST_DIR="${BUILD_DIR}/tests"
fi

echo

#----------------------------------------------------------------
echo "### Select Gauge group ###"
echo "1 [SU(3)] (default)"
echo "2 [SU(2)]"
echo "3 [General SU(N)]"
echo ""

echo -n "Gauge group? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	GAUGE_GROUP_FLAG="su3"
	echo "You have selected [SU(3)]"
	;;
    2)
	GAUGE_GROUP_FLAG="su2"
	echo "You have selected [SU(2)]"
	;;
    3)
	GAUGE_GROUP_FLAG="general"
	echo "You have selected [General SU(N)]"
	;;
    *)
	GAUGE_GROUP_FLAG="su3"
	echo "You have selected default value [SU(3)]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Select Target test name ###"
echo "1 [all] (default)"
echo "2 [Eigensolver]"
echo "3 [FFT]"
echo "4 [Gauge]"
echo "5 [GradientFlow]"
echo "6 [HMC]"
echo "7 [HotStart]"
echo "8 [IO]"
echo "9 [Mult]"
echo "10 [PolyakovLoop]"
echo "11 [QuarkNumberSusceptibility]"
echo "12 [RandomNumbers]"
echo "13 [Rational]"
echo "14 [SF_fAfP]"
echo "15 [Solver]"
echo "16 [Spectrum]"
echo "17 [TopologicalCharge]"
echo "18 [WilsonLoop]"
echo "19 [none]"
echo ""

echo -n "Test name? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	TEST_FLAG="all"
	echo "You have selected [all]"
	;;
    2)
	TEST_FLAG="Eigensolver"
	echo "You have selected [Eigensolver]"
	;;
    3)
	TEST_FLAG="FFT"
	echo "You have selected [FFT]"
	;;
    4)
	TEST_FLAG="Gauge"
	echo "You have selected [Gauge]"
	;;
    5)
	TEST_FLAG="GradientFlow"
	echo "You have selected [GradientFlow]"
	;;
    6)
	TEST_FLAG="HMC"
	echo "You have selected [HMC]"
	;;
    7)
	TEST_FLAG="HotStart"
	echo "You have selected [HotStart]"
	;;
    8)
	TEST_FLAG="IO"
	echo "You have selected [IO]"
	;;
    9)
	TEST_FLAG="Mult"
	echo "You have selected [Mult]"
	;;
    10)
	TEST_FLAG="PolyakovLoop"
	echo "You have selected [PolyakovLoop]"
	;;
    11)
	TEST_FLAG="QuarkNumberSusceptibility"
	echo "You have selected [QuarkNumberSusceptibility]"
	;;
    12)
	TEST_FLAG="RandomNumbers"
	echo "You have selected [RandomNumbers]"
	;;
    13)
	TEST_FLAG="Rational"
	echo "You have selected [Rational]"
	;;
    14)
	TEST_FLAG="SF_fAfP"
	echo "You have selected [SF_fAfP]"
	;;
    15)
	TEST_FLAG="Solver"
	echo "You have selected [Solver]"
	;;
    16)
	TEST_FLAG="Spectrum"
	echo "You have selected [Spectrum]"
	;;
    17)
	TEST_FLAG="TopologicalCharge"
	echo "You have selected [TopologicalCharge]"
	;;
    18)
	TEST_FLAG="WilsonLoop"
	echo "You have selected [WilsonLoop]"
	;;
    19)
	TEST_FLAG="none"
	echo "You have selected [none]"
	;;
    *)
	TEST_FLAG="all"
	echo "You have selected default value [all]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Test Manager ###"
echo "1 [Yes] (default)"
echo "2 [No]"
echo ""

echo -n "Use Test Manager? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	TESTMANAGER_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	TESTMANAGER_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	TESTMANAGER_FLAG="yes"
	echo "You have selected default value [Yes]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### External Library Options ###"
echo ""

#----------------------------------------------------------------
echo "### Extra library directory ###"
echo
echo "which directory to locate extra libraries?"
echo -n "[default: extra] > "
read VALUE

if [ -n "${VALUE}" ]; then
    EXTRA_DIR=${VALUE}
else
    EXTRA_DIR="extra"
fi

echo

#----------------------------------------------------------------
echo "### Use LIME library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use LIME library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_LIME_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_LIME_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_LIME_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use SFMT(SIMD Fast Mersenne Twister) library Option ###"
echo "## set [No] on FX1000"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use SFMT library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_SFMT_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_SFMT_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_SFMT_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use NTL library Option ###"
echo "## set [No] on FX1000"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use NTL library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_NTL_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_NTL_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_NTL_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use FFTW library Option ###"
echo "## set [No] on FX1000"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use FFTW library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_FFTW_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_FFTW_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_FFTW_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use YAML-CPP library Option ###"
echo "## set [No] on FX1000"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use YAML-CPP library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_YAMLCPP_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_YAMLCPP_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_YAMLCPP_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use tinyxml2 library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use tinyxml2 library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_TINYXML2_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_TINYXML2_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_TINYXML2_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use QWS library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use QWS library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_QWS_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_QWS_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_QWS_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use Alternative QXS library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Alternative QXS library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_ALT_QXS_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_ALT_QXS_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_ALT_QXS_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use QXS Architecture Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use QXS Architecture? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_QXS_ARCH_FLAG="general"
	echo "You have selected [general]"
	;;
    2)
	USE_QXS_ARCH_FLAG="acle"
	echo "You have selected [acle]"
	;;
    *)
	USE_QXS_ARCH_FLAG="general"
	echo "You have selected default value [general]"
	;;
esac

echo ""

#----------------------------------------------------------------

sed \
    -e "s/@@TARGET@@/$TARGET_FLAG/" \
    -e "s/@@USE_MPI@@/$MPI_FLAG/" \
    -e "s/@@USE_THREAD@@/$MULTITHREAD_FLAG/" \
    -e "s/@@USE_OPT_CODE@@/$OPT_CODE_FLAG/" \
    -e "s/@@USE_OPT_LEVEL@@/$OPT_LEVEL_FLAG/" \
    -e "s/@@USE_STATIC@@/$STATIC_FLAG/" \
    -e "s/@@USE_DEBUG@@/$DEBUG_FLAG/" \
    -e "s/@@USE_COMPLEX@@/$COMPLEX_FLAG/" \
    -e "s/@@USE_CPP11@@/$CPP11_FLAG/" \
    -e "s/@@USE_PYTHON@@/$PYTHON_FLAG/" \
    -e "s/@@USE_ALT_CODE@@/$ALT_CODE_FLAG/" \
    -e "s/@@USE_GAUGE_GROUP@@/$GAUGE_GROUP_FLAG/" \
    -e "s/@@TEST@@/$TEST_FLAG/" \
    -e "s/@@USE_TESTMANAGER@@/$TESTMANAGER_FLAG/" \
    -e "s/@@USE_LIME_LIBRARY@@/$USE_LIME_FLAG/" \
    -e "s/@@USE_SFMT_LIBRARY@@/$USE_SFMT_FLAG/" \
    -e "s/@@USE_NTL_LIBRARY@@/$USE_NTL_FLAG/" \
    -e "s/@@USE_FFTW_LIBRARY@@/$USE_FFTW_FLAG/" \
    -e "s/@@USE_YAMLCPP_LIBRARY@@/$USE_YAMLCPP_FLAG/" \
    -e "s/@@USE_TINYXML2_LIBRARY@@/$USE_TINYXML2_FLAG/" \
    -e "s/@@USE_QWS_LIBRARY@@/$USE_QWS_FLAG/" \
    -e "s/@@USE_ALT_QWS@@/$USE_ALT_QWS/" \
    -e "s/@@USE_QWS_ARCH@@/$USE_QWS_ARCH/" \
    -e "s,@@BUILD_DIR@@,$BUILD_DIR," \
    -e "s,@@TEST_DIR@@,$TEST_DIR," \
    -e "s,@@EXTRA_DIR@@,$EXTRA_DIR," \
Makefile.in > Makefile

#----------------------------------------------------------------
echo "Makefile is generated."

exit 0

