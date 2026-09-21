g++ -std=c++17 -O3 -march=native -flto -pthread -Ikmc_api kmc_merge_final.cpp ./kmc_api/*.cpp -o kmc_merge
mv kmc_merge ../
g++ -O2 -std=c++17 -pthread -o fisher_mt fisher_mt.cpp
mv fisher_mt ../
