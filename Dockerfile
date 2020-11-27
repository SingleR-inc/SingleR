FROM bioconductor/bioconductor_docker:devel

MAINTAINER infinite.monkeys.on.keyboards@gmail.com
LABEL authors="infinite.monkeys.on.keyboards@@gmail.com" \
    description="Docker image containing the SingleR package in a Bioconductor-devel container."

WORKDIR /home/build/package

COPY . /home/build/package 

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN Rscript -e "devtools::install('.', dependencies=TRUE, repos = BiocManager::repositories(), build_vignettes = TRUE)"
