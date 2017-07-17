#include <string>
#include "be_analyzer.h"

using namespace std;

int main (void) {
    const int R = 60, filted_n = 1, cleavage_length = 17, pri_len = 15; // 17bp + cleavage + 3bp + NGG for cas9
    //string grna_seq = "GAGTCCGAGCAGAAGAAGAA"; // 20bp for cas9
    //string wt_seq = "GGGCCTCCTGAGTTTCTCATCTGTGCCCCTCCCTCCCTGGCCCAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAACCACAAACCCACGAGGGCAGAGTGCTGCTTGCTGCTG"; // wt_seq for cas9
    //const char* fastq_joined_file = "/home/parkj/experiments/Be_Analyzer/data/7.fastqjoin"; // index.fastqjoin

    string grna_seq = "ATAGTCCAGGAGGCAGCCGA";
    string wt_seq = "CTGACGTGCCTCTCCCTCCCTCCAGGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATCTGCCTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGACTATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAGATCGCAAAGGTAATCAGGGAAGGGAGATACGGGGAGGGGAGATAAGG";
    const char* fastq_file1 = "/data/Sample_Cas_Analyzer_R1.fastq.gz"; // index.fastqjoin
    const char* fastq_file2 = "/data/Sample_Cas_Analyzer_R2.fastq.gz"; // index.fastqjoin

    // string control_file = "/home/parkj/experiments/Be_Analyzer/data/15.fastqjoin";

    BeAnalyzer ba = BeAnalyzer(R, filted_n, cleavage_length, pri_len, grna_seq, wt_seq);

    //ba.read_fastq_files(fastq_joined_file);
    ba.read_fastq_files(fastq_file1, fastq_file2);

    //ba.write_emboss();
    ba.run_alignment();
    ba.data_analyze();

    //ba.write_align();
    ba.write_count();
    //ba.write_substitution();
    puts("Hello World!");
    return 0;
}