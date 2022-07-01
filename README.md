
# Multivariate phenotype analysis enable genome-wide inference of mammalian gene function




## Installation

### Install Git and Docker
Follow instructions to install [Git](https://github.com/git-guides/install-git/) and [Docker](https://docs.docker.com/get-docker/) for your particular operating system.

### Clone Git repository

Choose the local directory where you want the repository cloned (let's call this "/path/to/your/git/repositories/"), and change directory to it. 


Linux:
```
MV_HOME="/path/to/your/local/git/repositories/"
cd $MV_HOME
```
Windows:
```
set MV_HOME="/path/to/your/local/git/repositories/"
cd %MV_HOME%
```
[Clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
our repository onto your local machine and checkout the paper release v2.0
```
git clone https://github.com/georgenicholson/multivariate_phenotype_data_and_code.git
cd multivariate_phenotype_data_and_code
git checkout tags/v2.0 -b my_local_branch
```
### Download and unzip results

Linux:
```
cd $MV_HOME/multivariate_phenotype_data_and_code
curl https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/download/v2.0/output.zip
unzip output.zip
```
Windows:
```
cd %MV_HOME%/multivariate_phenotype_data_and_code
curl https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/download/v2.0/output.zip
unzip output.zip
```
Alternatively, you can download output.zip from [here](https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/download/v2.0/output.zip), move the file to "%MV_HOME%/multivariate_phenotype_data_and_code", and unzip the file there, so that the results files are in "%MV_HOME%/multivariate_phenotype_data_and_code/output/global_results/".

### Pull or build Docker image

Pull the Docker image from [dockerhub](https://hub.docker.com/) which has been built from the Dockerfile in the root folder of the Git repository:
```
docker pull georgenicholson/multivariate_phenotype_data_and_code:v2.0
```
If instead you would like to build the Docker image locally, you can do this by running the following command:

Linux
```
docker build -t georgenicholson/multivariate_phenotype_data_and_code:v2.0 - < $MV_HOME/multivariate_phenotype_data_and_code/Dockerfile
```
Windows
```
docker build -t georgenicholson/multivariate_phenotype_data_and_code:v2.0 - < %MV_HOME%/multivariate_phenotype_data_and_code/Dockerfile
```



### Check installed location of data and results

The raw data were downloaded as part of the Git repository clone above, and should now be in [data/Data_all.RDS](data/Data_all.RDS).

The results files underlying the paper's Figures should be in subfolder "/output/global_results/"" once you have unzipped [output.zip](https://github.com/georgenicholson/multivariate_phenotype_data_and_code/releases/download/v2.01/output.zip) as described above.



## Generating paper figures

The script below, [scripts/05_generate_results.R](scripts/05_generate_results.R),  takes as input the raw data stored in [data/Data_all.RDS](data/Data_all.RDS) as well as the outputs of model fits, stored in "output/global_results/". Upon running the script below, the figures should be generated in the "figures/" folder.

Linux:
```
docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/05_generate_results.R
```
Windows:
```
docker run --entrypoint Rscript --rm --workdir /home -v %MV_HOME%/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/05_generate_results.R
```



## Demo run of software

The script below is a demo run on a data subset[^5], which takes approximately four minutes to complete and outputs results in "/output/methods_comparison".

[^5]:Demo run has: Method = ComposeMV, Data = impc, N = 100, P = 10, S = 2, K = 5

Linux:
```
docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R demo 1 1
```
Windows:
```
docker run --entrypoint Rscript --rm --workdir /home -v %MV_HOME%/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R demo 1 1
```

## Reproducing paper results

We list the analyses performed in the paper in Table A below. Each row of the table corresponds to an analysis based on a (Data, Method, Fold) combination, where "Fold" refers to cross-validation fold. Our method is labelled "ComposeMV" and we benchmark it alongside "mash" ([Urbut et al.](https://www.nature.com/articles/s41588-018-0268-8)), and "XD" ([Bovy et al.](https://www.jstor.org/stable/23024867?seq=1)), with our method building on both of these existing approaches. 

The analyses shown in Table A comprise the main analysis of our paper (row 5), our model checking analyses (rows 12-13), while the other rows comprise the benchmarking analysis presented in Tables 8-11 of our paper. N and P are a number of samples and measurement dimension of each data set. S and K are defined in the paper (S is the number of covariance matrices in the mixture model, and K is the dimensionality of the factor model). 


Each analysis performed on each of several cross validation folds (see # folds). Adequate memory allocations and approximate run times (single thread of an [Intel Xeon E5-1650 @ 3.20GHz](https://www.cpubenchmark.net/cpu.php?cpu=Intel+Xeon+E5-1650+%40+3.20GHz&id=1211)) are shown in Table A. 

The script to perform model fitting is "/scripts/01_model_fitting_wrapper.R" to ehich we provide three command line arguments:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(argument 1): "benchmark" to specify benchmarking analysis 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(argument 2): the analysis row of Table A as an integer 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(argument 3): the fold number as an integer 

### Reproducing one fold of the main analysis

The following command run the first fold of the main analysis of the paper (row 5 of Table A). The command-line arguments are "benchmark 5 1" and the analysis is run within the Docker image. Analysis outputs appear in "/output/methods_comparison/", with the analysis below taking approximately 1.8 hours to run.
 

Linux:
```
docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R benchmark 5 1
```
Windows:
```
docker run --entrypoint Rscript --rm --workdir /home -v %MV_HOME%/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R benchmark 5 1
```


### Reproducing all analyses

Below is pseudocode for reproducing all analyses, with "row" referring to a row of Table A and "fold" a cross-validation fold of the data set. Each (row, fold) analysis combination can be computed in parallel, so you can adapt the pseudocode to submit analyses as separate jobs to your server's workload management system. The total CPU time across all analyses is approximately four months, and the outputs are all generated in "/output/methods_comparison/".[^1] 
```
#!/bin/bash
for row in {1..24}
do
  for fold in {1..10} 
  do
    docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/01_model_fitting_wrapper.R benchmark $row $fold
  done
done
```
[^1]: Note from row 5 of Table A that the main analysis is run for 50 folds.




### Postprocessing analysis outputs into paper results

Once a full set of raw analysis outputs has been generated in "/output/methods_comparison/" the following four scripts are run sequentially to generate the results presented in the paper:

(1) Collect outputs together with "/scripts/02_collect_results.R"

(2) Fit factor model with "/scripts/03_estimate_global_factors.R"

(3) Perform hit calling "/scripts/04_calculate_hit_rates.R" 

(4) Generate figures and tables with "/scripts/05_generate_results.R" 

Each of the above scripts should be run within the Docker image, for example

Linux:
```
docker run --entrypoint Rscript --rm --workdir /home -v $MV_HOME/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/02_collect_results.R
```
Windows:
```
docker run --entrypoint Rscript --rm --workdir /home -v %MV_HOME%/multivariate_phenotype_data_and_code:/home georgenicholson/multivariate_phenotype_data_and_code:v2.0 --no-restore --no-save scripts/02_collect_results.R
```




### Table A: Analyses in paper


| Row | Method          | Data | N    | P   | S   | K  | # folds | RAM (GB) | CPU hrs / fold |
| :-- | :-------------- | :--- | :--- | :-- | :-- | :- | :------ | :------- | :------------- |
| 1   | XD              | impc | 2000 | 148 | 1   | NA | 10      | 2.3      | 3.1            |
| 2   | XD              | impc | 2000 | 148 | 2   | NA | 10      | 2.3      | 12.4           |
| 3   | mash            | impc | 2000 | 148 | 158 | NA | 10      | 6        | 11.2           |
| 4   | ComposeMV       | impc | 2000 | 148 | 1   | 15 | 10      | 2        | 1.7            |
| 5   | ComposeMV       | impc | 2000 | 148 | 1   | 20 | 50      | 2        | 1.8            |
| 6   | ComposeMV       | impc | 2000 | 148 | 1   | 30 | 10      | 2        | 1.8            |
| 7   | ComposeMV       | impc | 2000 | 148 | 1   | 40 | 10      | 2        | 1.8            |
| 8   | ComposeMV       | impc | 2000 | 148 | 2   | 15 | 10      | 2        | 25.7           |
| 9   | ComposeMV       | impc | 2000 | 148 | 2   | 20 | 10      | 2        | 28.4           |
| 10  | ComposeMV       | impc | 2000 | 148 | 2   | 30 | 10      | 2        | 26.8           |
| 11  | ComposeMV       | impc | 2000 | 148 | 2   | 40 | 10      | 2        | 20.9           |
| 12  | ComposeMV_N_500 | impc | 500  | 148 | 1   | 20 | 10      | 2.3      | 0.1            |
| 13  | ComposeMV_rand  | impc | 2000 | 148 | 1   | 20 | 10      | 2.3      | 0.3            |
| 14  | XD              | eqtl | 5000 | 44  | 1   | NA | 10      | 2.3      | 1.4            |
| 15  | XD              | eqtl | 5000 | 44  | 2   | NA | 10      | 2.3      | 3.7            |
| 16  | mash            | eqtl | 5000 | 44  | 54  | NA | 10      | 4.6      | 6.8            |
| 17  | ComposeMV       | eqtl | 5000 | 44  | 1   | 15 | 10      | 2        | 12.9           |
| 18  | ComposeMV       | eqtl | 5000 | 44  | 1   | 20 | 10      | 2        | 10.1           |
| 19  | ComposeMV       | eqtl | 5000 | 44  | 1   | 30 | 10      | 2        | 11.5           |
| 20  | ComposeMV       | eqtl | 5000 | 44  | 1   | 40 | 10      | 2        | 11.5           |
| 21  | ComposeMV       | eqtl | 5000 | 44  | 2   | 15 | 10      | 2        | 20.9           |
| 22  | ComposeMV       | eqtl | 5000 | 44  | 2   | 20 | 10      | 2        | 20.9           |
| 23  | ComposeMV       | eqtl | 5000 | 44  | 2   | 30 | 10      | 2        | 19.5           |
| 24  | ComposeMV       | eqtl | 5000 | 44  | 2   | 40 | 10      | 2        | 18             |
|     |

