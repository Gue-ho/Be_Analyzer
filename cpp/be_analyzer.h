//
// Created by parkj on 7/14/17.
//

#ifndef BE_ANALYZER_BE_ANALYZER_H
#define BE_ANALYZER_BE_ANALYZER_H

#include <string>
#include <map>
#include <vector>

using namespace std;

enum BeType {
    WT = 0,
    Ins = 1,
    Del = 2,
    Sub = 3,
    Other = 4
};

typedef map<string, int> seq_count;

struct sub_t {
    int cnt_A;
    int cnt_T;
    int cnt_G;
    int cnt_C;
    int& operator[](char sub) {
        if (sub == 'A') return cnt_A;
        else if (sub == 'T') return cnt_T;
        else if (sub == 'G') return cnt_G;
        else if (sub == 'C') return cnt_C;
    }
    const int& operator[](char sub) const {
        if (sub == 'A') return cnt_A;
        else if (sub == 'T') return cnt_T;
        else if (sub == 'G') return cnt_G;
        else if (sub == 'C') return cnt_C;
    }
};

class BeAnalyzer {
public:
    ~BeAnalyzer();
    BeAnalyzer(const int R, const int filted_n, const int cleavage_length, const int pri_len, string grna_seq, string wt_seq);
    void read_fastq_files(const char* fastq_file1, const char* fastq_file2);
    void read_fastq_files(const char* fastq_joined_file);
    int be_fastq_join(const char*, const char*);
    void run_alignment();
    //void write_emboss();
    void data_analyze();
    void write_align();
    void write_count();
    void write_substitution();
    void update_seq(string);

    int m_cleavage_position;
    int m_dir, m_cnt_all, m_cnt_pri, m_cnt_filt, m_cnt_sub, m_cnt_c_to_d;
    string m_wt_seq_sliced, m_pri_for, m_pri_back;
    int m_cnt_insertion, m_cnt_deletion, m_cnt_wt, m_cnt_others, m_R;

    vector<string> m_emboss_wt_list, m_emboss_seq_list, m_emboss_sym_list, m_sorted_list;

    seq_count m_seq_count;
    sub_t* m_pattern;

    int *m_cnt_position;
    int *m_cnt_c_to_d_position;

    int m_pri_len;

    vector<int> m_type;
};

#endif //BE_ANALYZER_BE_ANALYZER_H
