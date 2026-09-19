g++ -std=c++17 -O3 -march=native -flto -g -pthread kmc_merge_mt.cpp -o kmc_merge_mt -I ./kmc_api ./libkmc_core.a
