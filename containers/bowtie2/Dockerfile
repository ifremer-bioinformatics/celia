FROM alpine AS alpine-build-samtools

ARG SAMTOOLSVER=1.12
ARG BUILD_NCPU=1

RUN apk update; \
    apk upgrade; \
    apk --no-cache add \
    build-base \
    ncurses-dev ncurses-static \
    libbz2 \
    xz-dev \
    zlib-dev zlib-static \
    wget \
    bzip2-dev bzip2-static

WORKDIR /opt/src
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLSVER}/samtools-${SAMTOOLSVER}.tar.bz2
RUN tar -xjf samtools-${SAMTOOLSVER}.tar.bz2
WORKDIR /opt/src/samtools-${SAMTOOLSVER}
RUN ./configure
RUN make -j${BUILD_NCPU}
RUN strip samtools

FROM alpine AS alpine-build-bowtie2

ARG BOWTIE_VERSION=2.4.4
ARG BUILD_NCPU=2

RUN apk update; \
    apk add build-base wget zlib-dev zlib-static

WORKDIR /opt/src
RUN wget https://github.com/BenLangmead/bowtie2/archive/refs/tags/v${BOWTIE_VERSION}.tar.gz
RUN tar -xvf v${BOWTIE_VERSION}.tar.gz
WORKDIR /opt/src/bowtie2-${BOWTIE_VERSION}
RUN make -j${BUILD_NCPU} NO_TBB=1 EXTRA_CFLAGS=-static LDFLAGS=-static
RUN strip *-l *-s

FROM alpine

ARG BOWTIE_VERSION=2.4.4
ARG SAMTOOLS_VERSION=1.12

WORKDIR /
RUN apk add --no-cache bash
RUN apk add --no-cache perl
RUN apk add --no-cache python3
RUN apk add --no-cache libbz2 zlib libcurl xz-libs
COPY --from=alpine-build-bowtie2 /opt/src/bowtie2-${BOWTIE_VERSION}/bowtie2-* /usr/local/bin/
COPY --from=alpine-build-bowtie2 /opt/src/bowtie2-${BOWTIE_VERSION}/bowtie2 /usr/local/bin/
COPY --from=alpine-build-samtools /opt/src/samtools-${SAMTOOLS_VERSION}/samtools /usr/local/bin/samtools
