//
// Created by parkj on 7/18/17.
//
#include <emscripten.h>

#include <string>
#include <vector>
#include "be_analyzer.h"

extern "C" {
    void report_result(BeAnalyzer &ba) {
        int i;
        vector<string> types = {"WT", "Ins", "Del", "Sub", "Other"};
        for (i = 0; i < ba.m_sorted_list.size(); i++) {
            EM_ASM_({
                rtnval['align'].push([Pointer_stringify($0), Pointer_stringify($1), Pointer_stringify($2), Pointer_stringify($3), $4]);
            }, ba.m_emboss_wt_list[i].c_str(), ba.m_emboss_sym_list[i].c_str(), ba.m_emboss_seq_list[i].c_str(), types[ba.m_type[i]].c_str(), ba.m_seq_count[ba.m_sorted_list[i]]);
        }
    }

    void run_be_analyzer_paired(char *fastq_file1, char* fastq_file2, char* fullseq, char* rgenseq, int rev, int lrval, int nval) {
        int R = lrval, filted_n = nval, pri_len = 15, cleavage_length;
        string grna_seq(rgenseq);
        string wt_seq(fullseq);

        if (rev)
            cleavage_length = 4;
        else
            cleavage_length = grna_seq.size() - 3;

        BeAnalyzer ba = BeAnalyzer(R, filted_n, cleavage_length, pri_len, grna_seq, wt_seq);
        ba.read_fastq_files((string("/data/") + string(fastq_file1)).c_str(), (string("/data/") + string(fastq_file2)).c_str());
        ba.run_alignment();
        ba.data_analyze();
        ba.write_count();
        report_result(ba);
    }
    void run_be_analyzer_single(char *fastq_file, char* fullseq, char* rgenseq, int rev, int lrval, int nval) {
        int R = lrval, filted_n = nval, pri_len = 15, cleavage_length;
        string grna_seq(rgenseq);
        string wt_seq(fullseq);

        if (rev)
            cleavage_length = 4;
        else
            cleavage_length = grna_seq.size() - 3;

        BeAnalyzer ba = BeAnalyzer(R, filted_n, cleavage_length, pri_len, grna_seq, wt_seq);
        ba.read_fastq_files((string("/data/") + string(fastq_file)).c_str());
        ba.run_alignment();
        ba.data_analyze();
        report_result(ba);
    }

    void report_progress(float val, float total) {
        EM_ASM_({
            postMessage([0, $0]);
        }, val/total*100);
    }
};