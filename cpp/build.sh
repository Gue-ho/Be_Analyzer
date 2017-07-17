em++ -O3 -s USE_ZLIB=1 --std=c++11 main_v2.cpp be_fastq_join.cpp be_fastq-lib.cpp be_analyzer.cpp -o be_analyzer.js -I . -I fastq-join --pre-js needle.js --preload-file data
