FROM rocker/r-ver:4.1

RUN apt-get update

RUN apt-get install -y git

ADD https://raw.githubusercontent.com/georgenicholson/multivariate_phenotype_data_and_code/main/renv.lock ./

RUN R -e 'install.packages("renv")'

RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libpng-dev
RUN apt-get install -y libgsl-dev

RUN R -e 'renv::restore()'

COPY test_script.sh .
COPY generate_figures.sh .
COPY test_fits.R .
COPY demo.R .
COPY benchmark_1.sh .
COPY benchmark_2.sh .
COPY benchmark_3.sh .
COPY benchmark_1.R .
COPY benchmark_2.R .
COPY benchmark_3.R .
COPY run_demo.sh .

