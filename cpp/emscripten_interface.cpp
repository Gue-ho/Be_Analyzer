//
// Created by parkj on 7/18/17.
//
#include <emscripten.h>

#include <string>
#include <vector>
#include "be_analyzer.h"
#include "emscripten_progress.h"

extern "C" {
    void report_result(BeAnalyzer &ba) {
        int i, prev_prog = 0;
        vector<string> types = {"WT", "Ins", "Del", "Sub", "Other"};
        for (i = 0; i < ba.m_sorted_list.size(); i++) {
            int new_prog = int((i/(float)ba.m_sorted_list.size())*10);
            if (prev_prog < new_prog) {
                report_progress(i, ba.m_sorted_list.size(), 90, 10, "Preparing results...");
                prev_prog = new_prog;
            }
            EM_ASM_({
                rtnval['align'].push([Pointer_stringify($0), Pointer_stringify($1), Pointer_stringify($2), Pointer_stringify($3), $4]);
            }, ba.m_emboss_wt_list[i].c_str(), ba.m_emboss_sym_list[i].c_str(), ba.m_emboss_seq_list[i].c_str(), types[ba.m_type[i]].c_str(), ba.m_seq_count[ba.m_sorted_list[i]]);
        }
        EM_ASM_({rtnval['count'].push($0,$1,$2,$3,$4,$5,$6,$7,$8);},ba.m_cnt_all,ba.m_cnt_pri,ba.m_cnt_filt, ba.m_cnt_wt, ba.m_cnt_sub, ba.m_cnt_c_to_d, ba.m_cnt_insertion, ba.m_cnt_deletion, ba.m_cnt_others);
        EM_ASM_({rtnval['seq_info'].push(Pointer_stringify($0),Pointer_stringify($1),$2);},ba.m_wt_seq_sliced.c_str(),ba.m_grna_seq.c_str(), ba.m_R);
        EM_ASM({rtnval['sub'] = [[], [], [], []];});
        for (i = 0; i < ba.m_R*2; i++) {
            EM_ASM_({rtnval['sub'][0].push($0);},ba.m_pattern[i]['A']);
        }
        for (i = 0; i < ba.m_R*2; i++) {
            EM_ASM_({rtnval['sub'][1].push($0);},ba.m_pattern[i]['T']);
        }
        for (i = 0; i < ba.m_R*2; i++) {
            EM_ASM_({rtnval['sub'][2].push($0);},ba.m_pattern[i]['G']);
        }
        for (i = 0; i < ba.m_R*2; i++) {
            EM_ASM_({rtnval['sub'][3].push($0);},ba.m_pattern[i]['C']);
        }
        //report_progress(100, 100, 0, 100, "Done!");
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
};
