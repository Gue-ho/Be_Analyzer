em++ -O3 -s USE_ZLIB=1 --std=c++11 emscripten_interface.cpp be_fastq_join.cpp be_fastq-lib.cpp be_analyzer.cpp -o be-analyzer-core.js -I . -I fastq-join --pre-js needle.js -s EXPORTED_FUNCTIONS="['_run_be_analyzer_paired', '_run_be_analyzer_single']" --memory-init-file 0
cp be-analyzer-core.js* /home/baelab/Be_Analyzer/
