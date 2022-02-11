FROM jupyter/scipy-notebook:hub-1.4.2 AS build

LABEL maintainer="Nektar++ Development Team <nektar-users@imperial.ac.uk>"

USER root
RUN apt-get update && \
    apt-get install -y libgl1-mesa-dev xvfb openssh-client build-essential && \
    rm -rf /var/lib/apt/lists/*

USER $NB_UID
RUN conda install --quiet --yes jupyterlab=2.3.1 cmake=3.19.1 boost=1.74.0 occt=7.4.0 gmsh=4.6.0 pyvista=0.31.3 \
        arpack=3.7.0 fftw=3.3.8 scotch=6.0.9 openmpi=4.1.0 tinyxml=2.6.2 jupyter_contrib_nbextensions=0.5.1 \
        ptscotch=6.0.9 hdf5=1.10.6=mpi_openmpi* && \
    pip install --no-cache-dir nbgitpuller itkwidgets && \
    jupyter contrib nbextension install --sys-prefix && \
    jupyter nbextension enable --sys-prefix hide_input/main && \
    jupyter nbextension enable --sys-prefix hide_input_all/main && \
    jupyter nbextension enable --sys-prefix init_cell/main && \
    jupyter nbextension enable --sys-prefix equation-numbering/main && \
    jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib jupyterlab-datawidgets itkwidgets @jupyterlab/toc && \
    conda clean --all -f -y && \
    jupyter lab clean && \
    rm -rf "/home/${NB_USER}/.cache/yarn" && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

USER root
COPY . /tmp/nektar
RUN cd /tmp/nektar && mkdir build && cd build && \
    cmake -DNEKTAR_BUILD_SOLVERS=ON -DNEKTAR_SOLVER_DIFFUSION=OFF \
          -DNEKTAR_SOLVER_DUMMY=OFF -DNEKTAR_SOLVER_ELASTICITY=OFF \
          -DNEKTAR_SOLVER_MMF=OFF -DNEKTAR_SOLVER_CARDIAC_EP=OFF \
          -DNEKTAR_BUILD_UTILITIES=ON -DNEKTAR_BUILD_PYTHON=ON \
          -DNEKTAR_SOLVER_VORTEXWAVE=OFF -DNEKTAR_LIB_DIR=lib/x86_64-linux-gnu \
          -DNEKTAR_BUILD_TESTS=OFF -DNEKTAR_BUILD_UNIT_TESTS=OFF \
          -DNEKTAR_BUILD_DEMOS=OFF -DNEKTAR_USE_SCOTCH=ON -DNEKTAR_USE_MESHGEN=ON \
          -DNEKTAR_USE_ARPACK=ON -DNEKTAR_USE_FFTW=ON \
          -DNEKTAR_USE_MPI=ON -DNEKTAR_USE_HDF5=ON \
          -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release .. && \
    make -j$(nproc) install && \
    chown -R $NB_UID:$NB_UID /tmp/nektar && \
    cd /tmp/nektar/build && sudo -u $NB_USER make nekpy-install-system && cd && ldconfig && rm -R /tmp/nektar

USER root
COPY docker/nektar-workbook/entrypoint.sh /
RUN chmod +x /entrypoint.sh

USER $NB_USER
ENTRYPOINT ["/entrypoint.sh"]
