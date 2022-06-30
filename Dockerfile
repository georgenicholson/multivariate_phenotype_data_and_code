FROM bioconductor/bioconductor_docker:RELEASE_3_13
RUN apt-get update
RUN apt-get install -y git
ADD https://raw.githubusercontent.com/georgenicholson/multivariate_phenotype_data_and_code/main/renv.lock ./

ENV RENV_VERSION 0.15.1
ENV RENV_PATHS_LIBRARY /usr/local/lib/R/library
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libpng-dev
RUN apt-get install -y libgsl-dev

RUN R -e 'renv::activate()'
RUN R -e 'renv::restore()'
RUN R -e 'BiocManager::install(version = "3.13", ask = FALSE)'
RUN R -e 'BiocManager::install(pkgs = "AnnotationDbi", force = TRUE, version = "3.13")'
RUN R -e 'BiocManager::install(pkgs = "Mus.musculus", force = TRUE, version = "3.13")'
ENV R_PROFILE_USER="/"
ENV R_LIBS_USER="/usr/local/lib/R/library/R-4.1/x86_64-pc-linux-gnu"


