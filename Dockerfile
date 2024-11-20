# Build stage with Spack pre-installed and ready to be used
FROM spack/ubuntu-jammy:0.21 AS builder

# Install software from spack_devenv.yaml
RUN mkdir /opt/spack-environment
COPY spack.yaml /opt/spack-environment/spack.yaml
RUN cd /opt/spack-environment \
    && spack env activate . \
    && spack install --fail-fast \
    && spack gc -y

# Strip all the binaries
RUN find -L /opt/views/view/* -type f -exec readlink -f '{}' \; | \
    xargs file -i | \
    grep 'charset=binary' | \
    grep 'x-executable\|x-archive\|x-sharedlib' | \
    awk -F: '{print $1}' | xargs strip

# Modifications to the environment that are necessary to run
RUN cd /opt/spack-environment && \
    spack env activate --sh -d . > activate.sh

# Copy and compile domain_decomp
COPY . /decomp
RUN spack env activate /opt/spack-environment \
    && cd /decomp \
    && cmake -G "Unix Makefiles" -Bbuild -S. \
    && cmake --build build --config Release \
    && cp build/decomp /usr/bin/

# Bare OS image to run the installed executables
FROM ubuntu:22.04

# Copy necessary files from the builder stage
COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /decomp /decomp
COPY --from=builder /usr /usr
# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views
RUN ln -s /opt/views/view /opt/view

# Copy and set entrypoint
COPY entrypoint.sh /
RUN chmod a+x /entrypoint.sh
ENTRYPOINT [ "/entrypoint.sh" ]

# Set volume mount point for I/O
RUN mkdir /io
VOLUME /io
WORKDIR /io
