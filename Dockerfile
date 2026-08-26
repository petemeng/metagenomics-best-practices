FROM rocker/verse:4.4.1

ARG QUARTO_VERSION=1.9.38
ARG QUARTO_SHA256=bd0d73c1042fc5f719fb87753f0874c11344e516ee2cfbe4fe3ed4e4458c3988

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff-dev \
        libxml2-dev \
    && curl -fsSL \
        "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
        -o /tmp/quarto.deb \
    && echo "${QUARTO_SHA256}  /tmp/quarto.deb" | sha256sum -c - \
    && apt-get install -y --no-install-recommends /tmp/quarto.deb \
    && rm -f /tmp/quarto.deb \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project

COPY env/renv.lock env/renv.lock

RUN Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")' \
    && Rscript -e 'renv::restore(lockfile = "env/renv.lock", library = .libPaths()[1], prompt = FALSE)'

COPY . .

CMD ["quarto", "render"]
