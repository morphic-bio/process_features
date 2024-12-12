FROM debian:bookworm-slim as builder
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && apt-get -y install  build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libgomp1 libglib2.0-dev libcairo2-dev
ADD assignBarcodes.c /code/assignBarcodes.c
ADD Makefile /code/Makefile
ADD queue.h /code/queue.h
ADD queue.c /code/queue.c
ADD defines.h /code/defines.h
ADD plasma*.h  /code/.
WORKDIR /code
RUN make && mv assignBarcodes /usr/local/bin/assignBarcodes
RUN apt-get -y remove build-essential && apt-get -y autoremove

FROM debian:bookworm-slim
COPY --from=builder /usr/local/bin/assignBarcodes /usr/local/bin/assignBarcodes
#copy the zlib1g library and libglib2.0 libgomp.so.0
COPY --from=builder /usr/lib/x86_64-linux-gnu/* /usr/lib/x86_64-linux-gnu/.
# **Copy the font configuration files**
COPY --from=builder /etc/fonts /etc/fonts

# **Copy the system fonts**
COPY --from=builder /usr/share/fonts /usr/share/fonts