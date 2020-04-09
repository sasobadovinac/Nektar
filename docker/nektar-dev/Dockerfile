ARG ENV_IMAGE=nektarpp/nektar-env:latest
FROM $ENV_IMAGE AS build

LABEL maintainer="Nektar++ Development Team <nektar-users@imperial.ac.uk>"

ARG INSTALL_PREFIX=/usr/local
ARG NEKTAR_IMAGE=nektarpp/nektar:latest

USER root
COPY --from=${NEKTAR_IMAGE} ${INSTALL_PREFIX} ${INSTALL_PREFIX}

# Set up entrypoint for copying test files.
COPY docker/nektar/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh && \
    ln -s /usr/local/bin/docker-entrypoint.sh /
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

USER nektar
WORKDIR /home/nektar
CMD ["/bin/bash"]
