FROM debian:bookworm-slim AS builder
ARG VERSION=1.0.0
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && apt-get -y install  build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libgomp1 libglib2.0-dev curl libcairo2-dev libx11-dev libxext-dev libxrender-dev
RUN curl -L https://github.com/morphic-bio/process_features/archive/refs/tags/v${VERSION}.tar.gz | tar -xz
WORKDIR /process_features-${VERSION}
RUN make && mv assignBarcodes /usr/local/bin/assignBarcodes

FROM debian:bookworm-slim
COPY --from=builder /usr/local/bin/assignBarcodes /usr/local/bin/assignBarcodes
#copy the zlib1g library and libglib2.0 libgomp.so.0
#copy the libcairo2 libpixman-1 libfontconfig1 libfreetype6 libpng16 libxcb* lib-xrender  libraries
COPY --from=builder /usr/lib/x86_64-linux-gnu/libpng16.so.16 /usr/lib/x86_64-linux-gnu/libpng16.so.16
COPY --from=builder /usr/lib/x86_64-linux-gnu/libfontconfig.so.1 /usr/lib/x86_64-linux-gnu/libfontconfig.so.1
COPY --from=builder /usr/lib/x86_64-linux-gnu/libfreetype.so.6 /usr/lib/x86_64-linux-gnu/libfreetype.so.6
COPY --from=builder /usr/lib/x86_64-linux-gnu/libpixman-1.so.0 /usr/lib/x86_64-linux-gnu/libpixman-1.so.0
COPY --from=builder /usr/lib/x86_64-linux-gnu/libcairo.so.2 /usr/lib/x86_64-linux-gnu/libcairo.so.2
COPY --from=builder /usr/lib/x86_64-linux-gnu/libgomp.so.1 /usr/lib/x86_64-linux-gnu/libgomp.so.1
COPY --from=builder /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0 /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0
COPY --from=builder /usr/lib/x86_64-linux-gnu/libz.so.1 /usr/lib/x86_64-linux-gnu/libz.so.1
COPY --from=builder /usr/lib/x86_64-linux-gnu/libxcb* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libX* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libexpat.so.1 /usr/lib/x86_64-linux-gnu/libexpat.so.1
COPY --from=builder /usr/lib/x86_64-linux-gnu/libbrotlidec.so.1 /usr/lib/x86_64-linux-gnu/libbrotlidec.so.1
COPY --from=builder /usr/lib/x86_64-linux-gnu/libbrotlicommon.so.1 /usr/lib/x86_64-linux-gnu/libbrotlicommon.so.1
COPY --from=builder /usr/lib/x86_64-linux-gnu/libbsd.so.0 /usr/lib/x86_64-linux-gnu/libbsd.so.0

ENV OMP_NESTED=TRUE
ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
