#!/bin/bash

docker run --cpus=1 --cap-add SYS_ADMIN --cap-add DAC_READ_SEARCH -v /scratch:/scratch george_mv_test:latest /benchmark_1.sh > benchmark_1.sh.out.txt 2> benchmark_1.sh.err.txt &
docker run --cpus=3 --cap-add SYS_ADMIN --cap-add DAC_READ_SEARCH -v /scratch:/scratch george_mv_test:latest /benchmark_2.sh > benchmark_2.sh.out.txt 2> benchmark_2.sh.err.txt &
docker run --cpus=3 --cap-add SYS_ADMIN --cap-add DAC_READ_SEARCH -v /scratch:/scratch george_mv_test:latest /benchmark_3.sh > benchmark_3.sh.out.txt 2> benchmark_3.sh.err.txt &

wait;


