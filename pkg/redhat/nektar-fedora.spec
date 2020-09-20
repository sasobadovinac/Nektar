Name:           nektar++
Version:        %{_nektar_version}
Release:        %{_nektar_build_release}
Summary:        Nektar
License:        MIT
Group:          Development/Libraries/C and C++
Source0:        https://gitlab.nektar.info/nektar/nektar/-/archive/v%{version}/nektar-v%{version}.tar.gz
URL:            https://www.nektar.info/
Requires:       libnektar++ = %{version}
BuildRequires:  arpack-devel
BuildRequires:  blas-devel
BuildRequires:  boost-devel
BuildRequires:  boost-python3-devel
BuildRequires:  cmake
BuildRequires:  fftw-devel
BuildRequires:  gcc-c++
BuildRequires:  gcc-gfortran
BuildRequires:  lapack-devel
BuildRequires:  make
BuildRequires:  OCE-devel
BuildRequires:  petsc-devel
BuildRequires:  python3-devel
BuildRequires:  scotch-devel
BuildRequires:  tetgen-devel
BuildRequires:  tinyxml-devel
BuildRequires:  zlib-devel

%description
Nektar++ is a C++ framework for the development of computational solvers for
partial differential equations based on the spectral/hp element method.

This package contains the libraries which implement the spectral/hp element
discretisation and associated support classes.

#### System libraries
%package -n libnektar++
Summary:        Nektar++ spectral/hp element framework libraries
Group:          System/Libraries
Requires:       boost
Requires:       boost-python3
Requires:       tinyxml
Requires:       fftw
Requires:       arpack
Requires:       blas
Requires:       lapack
%description -n libnektar++
Nektar++ spectral/hp element framework libraries.

%package -n libnektar++-openmpi
BuildRequires:  openmpi
BuildRequires:  hdf5-openmpi-devel
BuildRequires:  petsc-openmpi-devel
BuildRequires:  ptscotch-openmpi-devel
Requires:       boost
Requires:       boost-python3
Requires:       tinyxml
Requires:       fftw
Requires:       arpack
Requires:       blas
Requires:       lapack
Requires:       hdf5-openmpi
Requires:       petsc-openmpi
Requires:       ptscotch-openmpi
Requires:       openmpi
Group:          System/Libraries
Summary:        OpenMPI variant of the Nektar++ spectral/hp element framework libraries
%description -n libnektar++-openmpi
OpenMPI variant of the Nektar++ spectral/hp element framework libraries

%package -n libnektar++-mpich
BuildRequires:  mpich
BuildRequires:  hdf5-mpich-devel
BuildRequires:  petsc-mpich-devel
BuildRequires:  ptscotch-mpich-devel
Requires:       boost
Requires:       boost-python3
Requires:       tinyxml
Requires:       fftw
Requires:       arpack
Requires:       blas
Requires:       lapack
Requires:       hdf5-mpich
Requires:       petsc-mpich
Requires:       ptscotch-mpich
Requires:       mpich
Group:          System/Libraries
Summary:        MPICH variant of the Nektar++ spectral/hp element framework libraries
%description -n libnektar++-mpich
MPICH variant of the Nektar++ spectral/hp element framework libraries

#### Development libraries
%package devel
Summary:        Development and header files for Nektar++
Group:          Development/Libraries/C and C++
Requires:       libnektar++ = %{version}
Requires:       arpack-devel
Requires:       blas-devel
Requires:       boost-devel
Requires:       boost-python3-devel
Requires:       cmake
Requires:       fftw-devel
Requires:       gcc-c++
Requires:       gcc-gfortran
Requires:       lapack-devel
Requires:       make
Requires:       OCE-devel
Requires:       petsc-devel
Requires:       python3-devel
Requires:       scotch-devel
Requires:       tetgen-devel
Requires:       tinyxml-devel
Requires:       zlib-devel
%description devel
Development and header files for Nektar++

%package openmpi-devel
Summary:        Development and header files for Nektar++ (OpenMPI variant)
Group:          Development/Libraries/C and C++
Requires:       libnektar++-openmpi = %{version}
Requires:       arpack-devel
Requires:       blas-devel
Requires:       boost-devel
Requires:       boost-python3-devel
Requires:       cmake
Requires:       fftw-devel
Requires:       gcc-c++
Requires:       gcc-gfortran
Requires:       lapack-devel
Requires:       make
Requires:       OCE-devel
Requires:       petsc-openmpi-devel
Requires:       hdf5-openmpi-devel
Requires:       python3-devel
Requires:       ptscotch-openmpi-devel
Requires:       tetgen-devel
Requires:       tinyxml-devel
Requires:       zlib-devel
%description openmpi-devel
Development and header files for Nektar++ (OpenMPI variant)

%package mpich-devel
Summary:        Development and header files for Nektar++ (MPICH variant)
Group:          Development/Libraries/C and C++
Requires:       libnektar++-mpich = %{version}
Requires:       arpack-devel
Requires:       blas-devel
Requires:       boost-devel
Requires:       boost-python3-devel
Requires:       cmake
Requires:       fftw-devel
Requires:       gcc-c++
Requires:       gcc-gfortran
Requires:       lapack-devel
Requires:       make
Requires:       OCE-devel
Requires:       petsc-mpich-devel
Requires:       hdf5-mpich-devel
Requires:       python3-devel
Requires:       ptscotch-mpich-devel
Requires:       tetgen-devel
Requires:       tinyxml-devel
Requires:       zlib-devel
%description mpich-devel
Development and header files for Nektar++ (MPICH variant)

#### Python bindings
%package python3
Summary:        Nektar++ Python 3 interface library
Requires:       libnektar++ = %{version}
%description python3
Nektar++ Python 3 interface library

%package python3-openmpi
Summary:        Nektar++ Python 3 interface library (OpenMPI variant)
Requires:       libnektar++-openmpi = %{version}
Requires:       python3-openmpi
%description python3-openmpi
Nektar++ Python 3 interface library (OpenMPI variant)

%package python3-mpich
Summary:        Nektar++ Python 3 interface library (MPICH variant)
Requires:       libnektar++ = %{version}
Requires:       python3-mpich
%description python3-mpich
Nektar++ Python 3 interface library (MPICH variant)

#### Documentation
%package doc
Summary:        Documentation for Nektar++
Group:          Documentation
BuildRequires:  doxygen
BuildRequires:  ImageMagick
BuildRequires:  graphviz
BuildRequires:  texlive texlive-import texlive-lstaddons texlive-bclogo
BuildRequires:  texlive-mdframed texlive-standalone
%description doc
Documentation for Nektar++

#### Utilities

## NekMesh (serial only, no MPI support)
%package -n nekmesh
Summary:        NekMesh high-order mesh generation framework
Group:          Scientific
Requires:       libnektar++ = %{version}
%description -n nekmesh
NekMesh high-order mesh generation framework

## FieldConvert
%package fieldconvert
Summary:        FieldConvert post-processing utility for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description fieldconvert
FieldConvert post-processing utility for Nektar++

%package openmpi-fieldconvert
Summary:        FieldConvert post-processing utility for Nektar++
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-fieldconvert
FieldConvert post-processing utility for Nektar++

%package mpich-fieldconvert
Summary:        FieldConvert post-processing utility for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description mpich-fieldconvert
FieldConvert post-processing utility for Nektar++

#### Solvers

## AcousticSolver
%package acoustic-solver
Summary:        Acoustic solver for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description acoustic-solver
Acoustic solver for Nektar++

%package openmpi-acoustic-solver
Summary:        Acoustic solver for Nektar++ (OpenMPI)
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-acoustic-solver
Acoustic solver for Nektar++ (OpenMPI)

%package mpich-acoustic-solver
Summary:        Acoustic solver for Nektar++ (MPICH)
Group:          Scientific
Requires:       libnektar++-mpich = %{version}
%description mpich-acoustic-solver
Acoustic solver for Nektar++ (MPICH)

## ADRSolver
%package adr-solver
Summary:        Advection-diffusion-reaction solver for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description adr-solver
Advection-diffusion-reaction solver for Nektar++

%package openmpi-adr-solver
Summary:        Advection-diffusion-reactionsolver for Nektar++ (OpenMPI)
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-adr-solver
Advection-diffusion-reactionsolver for Nektar++ (OpenMPI)

%package mpich-adr-solver
Summary:        Advection-diffusion-reactionsolver for Nektar++ (MPICH)
Group:          Scientific
Requires:       libnektar++-mpich = %{version}
%description mpich-adr-solver
Advection-diffusion-reactionsolver for Nektar++ (MPICH)

## CardiacEPSolver
%package cardiacep-solver
Summary:        Cardiac electrophysiology solver for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description cardiacep-solver
Cardiac electrophysiology solver for Nektar++

%package openmpi-cardiacep-solver
Summary:        Cardiac electrophysiology solver for Nektar++ (OpenMPI)
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-cardiacep-solver
Cardiac electrophysiology solver for Nektar++ (OpenMPI)

%package mpich-cardiacep-solver
Summary:        Cardiac electrophysiology solver for Nektar++ (MPICH)
Group:          Scientific
Requires:       libnektar++-mpich = %{version}
%description mpich-cardiacep-solver
Cardiac electrophysiology solver for Nektar++ (MPICH)

## CompressibleFlowSolver
%package compressibleflow-solver
Summary:        Compressible flow solver for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description compressibleflow-solver
Compressible flow solver for Nektar++

%package openmpi-compressibleflow-solver
Summary:        Compressible flow solver for Nektar++ (OpenMPI)
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-compressibleflow-solver
Compressible flow solver for Nektar++ (OpenMPI)

%package mpich-compressibleflow-solver
Summary:        Compressible flow solver for Nektar++ (MPICH)
Group:          Scientific
Requires:       libnektar++-mpich = %{version}
%description mpich-compressibleflow-solver
Compressible flow solver for Nektar++ (MPICH)

## IncNavierStokesSolver
%package incnavierstokes-solver
Summary:        Incompressible flow solver for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description incnavierstokes-solver
Incompressible flow solver for Nektar++

%package openmpi-incnavierstokes-solver
Summary:        Incompressible flow solver for Nektar++ (OpenMPI)
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-incnavierstokes-solver
Incompressible flow solver for Nektar++ (OpenMPI)

%package mpich-incnavierstokes-solver
Summary:        Incompressible flow solver for Nektar++ (MPICH)
Group:          Scientific
Requires:       libnektar++-mpich = %{version}
%description mpich-incnavierstokes-solver
Incompressible flow solver for Nektar++ (MPICH)

## PulseWaveSolver
%package pulsewave-solver
Summary:        Pulse wave solver for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description pulsewave-solver
Pulse wave flow solver for Nektar++

%package openmpi-pulsewave-solver
Summary:        Pulse wave flow solver for Nektar++ (OpenMPI)
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-pulsewave-solver
Pulse wave flow solver for Nektar++ (OpenMPI)

%package mpich-pulsewave-solver
Summary:        Pulse wave flow solver for Nektar++ (MPICH)
Group:          Scientific
Requires:       libnektar++-mpich = %{version}
%description mpich-pulsewave-solver
Pulse wave flow solver for Nektar++ (MPICH)

## ShallowWaterSolver
%package shallowwater-solver
Summary:        Shallow water solver for Nektar++
Group:          Scientific
Requires:       libnektar++ = %{version}
%description shallowwater-solver
Shallow water flow solver for Nektar++

%package openmpi-shallowwater-solver
Summary:        Shallow water flow solver for Nektar++ (OpenMPI)
Group:          Scientific
Requires:       libnektar++-openmpi = %{version}
%description openmpi-shallowwater-solver
Shallow water flow solver for Nektar++ (OpenMPI)

%package mpich-shallowwater-solver
Summary:        Shallow water flow solver for Nektar++ (MPICH)
Group:          Scientific
Requires:       libnektar++-mpich = %{version}
%description mpich-shallowwater-solver
Shallow water flow solver for Nektar++ (MPICH)

%prep
%setup -qn nektar-v%{version}

%define dobuild()                                           \
mkdir $MPI_COMPILER;                                        \
cd $MPI_COMPILER;                                           \
%{cmake}                                                  \\\
    -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX                \\\
    -DCMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES=/usr/include \\\
    -DNEKTAR_ERROR_ON_WARNINGS=OFF                        \\\
    -DNEKTAR_LIB_DIR=$NEKTAR_LIBDIR                       \\\
    -DNEKTAR_INCLUDE_ROOT=$NEKTAR_INCLUDE_ROOT            \\\
    -DNEKTAR_BUILD_DEMOS=OFF                              \\\
    -DNEKTAR_BUILD_TESTS=OFF                              \\\
    -DNEKTAR_BUILD_UNIT_TESTS=OFF                         \\\
    -DNEKTAR_BUILD_SOLVERS=ON                             \\\
    -DNEKTAR_SOLVER_DIFFUSION=OFF                         \\\
    -DNEKTAR_SOLVER_DUMMY=OFF                             \\\
    -DNEKTAR_SOLVER_DIFFUSION=OFF                         \\\
    -DNEKTAR_SOLVER_IMAGE_WARPING=OFF                     \\\
    -DNEKTAR_SOLVER_ELASTICITY=OFF                        \\\
    -DNEKTAR_SOLVER_DIFFUSION=OFF                         \\\
    -DNEKTAR_SOLVER_MMF=OFF                               \\\
    -DNEKTAR_SOLVER_VORTEXWAVE=OFF                        \\\
    -DNEKTAR_BUILD_PYTHON=ON                              \\\
    -DNEKTAR_USE_MESHGEN=ON                               \\\
    -DNEKTAR_USE_MPI=$MPI_ON                              \\\
    -DNEKTAR_USE_ARPACK=ON                                \\\
    -DNEKTAR_USE_FFTW=ON                                  \\\
    -DNEKTAR_USE_HDF5=$MPI_ON                             \\\
    -DNEKTAR_USE_CCM=OFF                                  \\\
    -DNEKTAR_USE_PYTHON3=ON                               \\\
    .. ;                                                    \
make %{?_smp_mflags} ; \
cd ..

%build

%undefine _hardened_build

# Build serial version, dummy arguments
MPI_COMPILER=serial MPI_SUFFIX= MPI_ON=OFF NEKTAR_LIBDIR=%{_lib} NEKTAR_INCLUDE_ROOT=%{_prefix}/include INSTALL_PREFIX=%{_prefix} %dobuild

# Build documentation.
cd serial && cmake -DNEKTAR_BUILD_DOC=ON .
make user-guide-pdf developer-guide-pdf doc
cd ..

# Make sure module paths are loaded
. /etc/profile.d/00-modulepath.sh;

# Build OpenMPI version
%{_openmpi_load}
MPI_ON=ON NEKTAR_LIBDIR=lib NEKTAR_INCLUDE_ROOT=%{_prefix}/include/$MPI_COMPILER INSTALL_PREFIX=%{_prefix}/lib64/openmpi %dobuild
%{_openmpi_unload}

# Build mpich version
%{_mpich_load}
MPI_ON=ON NEKTAR_LIBDIR=lib NEKTAR_INCLUDE_ROOT=%{_prefix}/include/$MPI_COMPILER INSTALL_PREFIX=%{_prefix}/lib64/mpich %dobuild
%{_mpich_unload}

%install
# Install serial version
make -C serial install DESTDIR=%{buildroot} INSTALL="install -p" CPPROG="cp -p"

# Install serial NekPy library
cd serial
%{__python3} setup.py install --root=%{buildroot} --install-purelib=%{python3_sitearch}
cd ..

# Make sure module paths are loaded, again.
. /etc/profile.d/00-modulepath.sh;

# Install OpenMPI version
%{_openmpi_load}
cd $MPI_COMPILER
make install DESTDIR=%{buildroot}
%{__python3} setup.py install --root=%{buildroot} --install-purelib=%{python3_sitearch}/openmpi
mv %{buildroot}/usr/lib64/openmpi/include %{buildroot}/usr/include/$MPI_COMPILER
cd ..
%{_openmpi_unload}

# Install MPICH version
%{_mpich_load}
cd $MPI_COMPILER
make install DESTDIR=%{buildroot}
%{__python3} setup.py install --root=%{buildroot} --install-purelib=%{python3_sitearch}/mpich
mv %{buildroot}/usr/lib64/mpich/include %{buildroot}/usr/include/$MPI_COMPILER
cd ..
%{_mpich_unload}

# Clean up temporary Python files
rm -rf %{buildroot}/root

# Remove MPI NekMesh executables
rm %{buildroot}/usr/lib64/mpich/bin/NekMesh
rm %{buildroot}/usr/lib64/openmpi/bin/NekMesh

#### Files for RPM packages

%files -n libnektar++
%defattr(-,root,root,-)
%{_libdir}/*.so.*

%files -n libnektar++-openmpi
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so.*

%files -n libnektar++-mpich
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so.*

%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
%{_libdir}/nektar++
%{_includedir}/*

%files openmpi-devel
%defattr(-,root,root,-)
%{_libdir}/openmpi/lib/*.so
%{_libdir}/openmpi/lib/nektar++
%{_includedir}/openmpi-x86_64/*

%files mpich-devel
%defattr(-,root,root,-)
%{_libdir}/mpich/lib/*.so
%{_libdir}/mpich/lib/nektar++
%{_includedir}/mpich-x86_64/*

%files python3
%{python3_sitearch}/NekPy
%{python3_sitearch}/NekPy*.egg-info

%files python3-openmpi
%{python3_sitearch}/openmpi/NekPy
%{python3_sitearch}/openmpi/NekPy*.egg-info

%files python3-mpich
%{python3_sitearch}/mpich/NekPy
%{python3_sitearch}/mpich/NekPy*.egg-info

## Documentation
%files doc
%{_docdir}/nektar++

## Utilities
%files -n nekmesh
%{_bindir}/NekMesh

%files fieldconvert
%{_bindir}/FieldConvert

%files openmpi-fieldconvert
%{_libdir}/openmpi/bin/FieldConvert

%files mpich-fieldconvert
%{_libdir}/mpich/bin/FieldConvert

## Solvers
%files acoustic-solver
%{_bindir}/AcousticSolver
%{_bindir}/APESolver
%files openmpi-acoustic-solver
%{_libdir}/openmpi/bin/AcousticSolver
%{_libdir}/openmpi/bin/APESolver
%files mpich-acoustic-solver
%{_libdir}/mpich/bin/AcousticSolver
%{_libdir}/mpich/bin/APESolver

%files adr-solver
%{_bindir}/ADRSolver
%files openmpi-adr-solver
%{_libdir}/openmpi/bin/ADRSolver
%files mpich-adr-solver
%{_libdir}/mpich/bin/ADRSolver

%files cardiacep-solver
%{_bindir}/CardiacEPSolver
%{_bindir}/PrePacing
%files openmpi-cardiacep-solver
%{_libdir}/openmpi/bin/CardiacEPSolver
%{_libdir}/openmpi/bin/PrePacing
%files mpich-cardiacep-solver
%{_libdir}/mpich/bin/CardiacEPSolver
%{_libdir}/mpich/bin/PrePacing

%files compressibleflow-solver
%{_bindir}/CompressibleFlowSolver
%{_bindir}/CompressibleBL
%{_bindir}/ExtractSurface2DCFS
%{_bindir}/ExtractSurface3DCFS
%files openmpi-compressibleflow-solver
%{_libdir}/openmpi/bin/CompressibleFlowSolver
%{_libdir}/openmpi/bin/CompressibleBL
%{_libdir}/openmpi/bin/ExtractSurface2DCFS
%{_libdir}/openmpi/bin/ExtractSurface3DCFS
%files mpich-compressibleflow-solver
%{_libdir}/mpich/bin/CompressibleFlowSolver
%{_libdir}/mpich/bin/CompressibleBL
%{_libdir}/mpich/bin/ExtractSurface2DCFS
%{_libdir}/mpich/bin/ExtractSurface3DCFS

%files incnavierstokes-solver
%{_bindir}/IncNavierStokesSolver
%{_bindir}/AddModeTo2DFld
%{_bindir}/Aliasing
%{_bindir}/CFLStep
%{_bindir}/ExtractMeanModeFromHomo1DFld
%{_bindir}/Fld2DTo2D5
%{_bindir}/FldAddFalknerSkanBL
%{_bindir}/NonLinearEnergy
%files openmpi-incnavierstokes-solver
%{_libdir}/openmpi/bin/IncNavierStokesSolver
%{_libdir}/openmpi/bin/AddModeTo2DFld
%{_libdir}/openmpi/bin/Aliasing
%{_libdir}/openmpi/bin/CFLStep
%{_libdir}/openmpi/bin/ExtractMeanModeFromHomo1DFld
%{_libdir}/openmpi/bin/Fld2DTo2D5
%{_libdir}/openmpi/bin/FldAddFalknerSkanBL
%{_libdir}/openmpi/bin/NonLinearEnergy
%files mpich-incnavierstokes-solver
%{_libdir}/mpich/bin/IncNavierStokesSolver
%{_libdir}/mpich/bin/AddModeTo2DFld
%{_libdir}/mpich/bin/Aliasing
%{_libdir}/mpich/bin/CFLStep
%{_libdir}/mpich/bin/ExtractMeanModeFromHomo1DFld
%{_libdir}/mpich/bin/Fld2DTo2D5
%{_libdir}/mpich/bin/FldAddFalknerSkanBL
%{_libdir}/mpich/bin/NonLinearEnergy

%files pulsewave-solver
%{_bindir}/PulseWaveSolver
%{_bindir}/Fld2Tecplot
%files openmpi-pulsewave-solver
%{_libdir}/openmpi/bin/PulseWaveSolver
%{_libdir}/openmpi/bin/Fld2Tecplot
%files mpich-pulsewave-solver
%{_libdir}/mpich/bin/PulseWaveSolver
%{_libdir}/mpich/bin/Fld2Tecplot

%files shallowwater-solver
%{_bindir}/ShallowWaterSolver
%files openmpi-shallowwater-solver
%{_libdir}/openmpi/bin/ShallowWaterSolver
%files mpich-shallowwater-solver
%{_libdir}/mpich/bin/ShallowWaterSolver

