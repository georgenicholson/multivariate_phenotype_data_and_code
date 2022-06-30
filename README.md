# Multivariate phenotype analysis enable genome-wide inference of mammalian gene function




## Installation

### Install Git and Docker
Follow instructions to install [Git](https://github.com/git-guides/install-git/) and [Docker](https://docs.docker.com/get-docker/) for your particular operating system.

### Clone the Git repository

Choose the local directory where you want the repository cloned (let's call this “/path/to/your/git/repos/”), and change directory to it. 


Linux:
```
MV_HOME="/path/to/your/local/git/repositories/"
cd $MV_HOME
```
Windows:
```
set MV_HOME=“/path/to/your/local/git/repositories/”
cd %MV_HOME%
```
[Clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
our repository onto your local machine
```
git clone https://github.com/georgenicholson/multivariate_phenotype_data_and_code.git

```
### Download and unzip the results files



Linux:
```
cd $MV_HOME/multivariate_phenotype_data_and_code
wget https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/download/v1.1/output.zip
unzip output.zip
```
Windows users, please download the results in [output.zip](https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/download/v1.1/output.zip), move the file to "%MV_HOME%/multivariate_phenotype_data_and_code", and unzip the file there.


### Pull the Docker image

Pull the Docker image from [dockerhub](https://hub.docker.com/) which has been created from the Dockerfile in the root folder of the git repository

```
docker pull georgenicholson/multivariate_phenotype_data_and_code:v2.0
```

### Data and results

The raw data were downloaded as part of the Git repository clone above, and should now be in [/data/Data_all.RDS](/data/Data_all.RDS).

The results files underlying the paper's Figures should be in subfolder [/output/global_results/](/output/global_results/) once you have unzipped [output.zip](https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/download/v1.1/output.zip) as described above.



## Generating the figures from the paper

The script below, [scripts/05_generate_results.R](scripts/05_generate_results.R),  takes as input the raw data stored in [/data/Data_all.RDS](data/Data_all.RDS) as well as the outputs of model fits, stored in [/output/global_results/](output/global_results/).

Linux:
```
docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/05_generate_results.R
```
Windows:
```
docker run --entrypoint Rscript --rm --workdir /home -v %MV_HOME%/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/05_generate_results.R
```



The Figures should now have been generated in the [/figures/](/figures/) folder.

## Demo run of software

Linux:
```
docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R demo 1 1
```
Windows:
```
docker run --entrypoint Rscript --rm --workdir /home -v %MV_HOME%/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R demo 1 1
```


## Reproducing our paper's results

We represent the analyses performed in the paper in Table A below. Each row of the table corresponds to an analysis based on a (Data, Method) pair. Our method is labelled "ComposeMV" and we benchmark it alongside "mash" ([Urbut et al.](https://www.nature.com/articles/s41588-018-0268-8)), and "XD" ([Bovy et al.](https://www.jstor.org/stable/23024867?seq=1))---our methods build on both of these existing approaches. 

The analyses shown in Table A comprise the main analysis of our paper (row 5), our model checking analyses (rows 12-13), while the other rows comprise the benchmarking analysis (Tables 8-11 of our paper). N and P are a number of samples and dimensionality of each data set. S and K are defined in the paper (S is the number of covariance matrices in the mixture, and K is the dimensionality of the factor model). 


Each analysis performed on each of several cross validation folds (see # folds). Adequate memory allocations and approximate run times (single thread of an [Intel Xeon E5-1650 @ 3.20GHz](https://www.cpubenchmark.net/cpu.php?cpu=Intel+Xeon+E5-1650+%40+3.20GHz&id=1211)) are shown in Table A. 

The script to perform model fitting is "/scripts/01_model_fitting_wrapper.R" and we provide three command line arguments: (i) "benchmark" to specify benchmarking analysis; (ii) the analysis row as an integer; (iii) the fold number as an integer. For example to run the first fold of the main analysis of the paper (row 5 of Table A), the command line arguments are "benchmark 5 1" and the analysis is run within the Docker image as follows:

Linux:
```
docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R benchmark 5 1
```
Windows:
```
docker run --entrypoint Rscript --rm --workdir /home -v %MV_HOME%/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R benchmark 5 1
```


### Table A: benchmarking analyses


| Row | Data | N    | P   | Method          | S   | K  | N folds | RAM (GB) | CPU hrs / fold |
| :-- | :--- | :--- | :-- | :-------------- | :-- | :- | :------ | :------- | :------------- |
| 1   | impc | 2000 | 148 | XD              | 1   | NA | 10      | 2.3      | 3.1            |
| 2   | impc | 2000 | 148 | XD              | 2   | NA | 10      | 2.3      | 12.4           |
| 3   | impc | 2000 | 148 | mash            | 158 | NA | 10      | 6        | 11.2           |
| 4   | impc | 2000 | 148 | ComposeMV       | 1   | 15 | 10      | 2        | 1.7            |
| 5   | impc | 2000 | 148 | ComposeMV       | 1   | 20 | 50      | 2        | 1.8            |
| 6   | impc | 2000 | 148 | ComposeMV       | 1   | 30 | 10      | 2        | 1.8            |
| 7   | impc | 2000 | 148 | ComposeMV       | 1   | 40 | 10      | 2        | 1.8            |
| 8   | impc | 2000 | 148 | ComposeMV       | 2   | 15 | 10      | 2        | 25.7           |
| 9   | impc | 2000 | 148 | ComposeMV       | 2   | 20 | 10      | 2        | 28.4           |
| 10  | impc | 2000 | 148 | ComposeMV       | 2   | 30 | 10      | 2        | 26.8           |
| 11  | impc | 2000 | 148 | ComposeMV       | 2   | 40 | 10      | 2        | 20.9           |
| 12  | impc | 500  | 148 | ComposeMV_N_500 | 1   | 20 | 10      | 2.3      | 0.1            |
| 13  | impc | 2000 | 148 | ComposeMV_rand  | 1   | 20 | 10      | 2.3      | 0.3            |
| 14  | eqtl | 5000 | 44  | XD              | 1   | NA | 10      | 2.3      | 1.4            |
| 15  | eqtl | 5000 | 44  | XD              | 2   | NA | 10      | 2.3      | 3.7            |
| 16  | eqtl | 5000 | 44  | mash            | 54  | NA | 10      | 4.6      | 6.8            |
| 17  | eqtl | 5000 | 44  | ComposeMV       | 1   | 15 | 10      | 2        | 12.9           |
| 18  | eqtl | 5000 | 44  | ComposeMV       | 1   | 20 | 10      | 2        | 10.1           |
| 19  | eqtl | 5000 | 44  | ComposeMV       | 1   | 30 | 10      | 2        | 11.5           |
| 20  | eqtl | 5000 | 44  | ComposeMV       | 1   | 40 | 10      | 2        | 11.5           |
| 21  | eqtl | 5000 | 44  | ComposeMV       | 2   | 15 | 10      | 2        | 20.9           |
| 22  | eqtl | 5000 | 44  | ComposeMV       | 2   | 20 | 10      | 2        | 20.9           |
| 23  | eqtl | 5000 | 44  | ComposeMV       | 2   | 30 | 10      | 2        | 19.5           |
| 24  | eqtl | 5000 | 44  | ComposeMV       | 2   | 40 | 10      | 2        | 18             |
|     |
