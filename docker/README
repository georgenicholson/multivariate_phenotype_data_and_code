# Docker Reproducability Instructions

This directory contains the files required to build a docker image of the code for the paper 'Multivariate phenotype analysis enables genome-wide inference of mammalian gene function'.

## Run the image from docker hub

NOTE: This currently requires access to the MRC Harwell docker hub (cutter).  This shall be pushed to a public docker hub.

### Install docker
Follow instructions at https://docs.docker.com/get-docker/ for your operating system

### Clone the repository
To get started, choose the local directory where you want the repository (let's call this <your_path>), and [clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
this repository onto your local machine in folder "<your_path>/multivariate_phenotype_data_and_code"
```
cd <your_path>
git clone https://github.com/georgenicholson/multivariate_phenotype_data_and_code.git
```

### Run demo
```
docker run --name hughtest --rm -v <your_path>/multivariate_phenotype_data_and_code:/multivariate_phenotype_data_and_code cutter:5000/george_mv_test /run_demo.sh
```

### Run benchmarks
```
nohup ./run_benchmarks_14_12_23.sh
```

## Build your own image
To build your own image locally:
```
cd <your_path>
docker build -t george_mv_test .
```

