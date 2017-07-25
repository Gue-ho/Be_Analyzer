//
// Created by parkj on 7/17/17.
//

#include <string>
#include <algorithm>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream> //sort
#include <set>
#include <functional>
#include "be_fastq-lib.h"
#include "be_fastq-join.h"
#include "be_analyzer.h"
#include "emscripten_progress.h"

using namespace std;

struct less_second {
    typedef pair<string, int> type;
    bool operator () (type const& a, type const &b) const {
        return a.second > b.second;
    }
};

vector<string> sorted_filtered_map(seq_count m_out_seq, int crit, int *filtered_cnt) {
    vector<pair<string, int> > cp;

    *filtered_cnt = 0;
    for (seq_count::iterator iter=m_out_seq.begin(); iter != m_out_seq.end(); iter++) {
        if (iter->second > crit) {
            cp.push_back(*iter);
            *filtered_cnt += iter->second;
        }
    }

    sort(cp.begin(), cp.end(), less_second());

    vector<string> rtn;
    for (int i=0; i<cp.size(); i++)
        rtn.push_back(cp[i].first);

    return rtn;
}

string rev_comp(string &a) {
    //string a must be uppercase letters
    int i;
    string s;

    for (i = 0; i<a.size(); i++) {
        if (a[i] == 'A') a[i] = 'T';
        else if (a[i] == 'T') a[i] = 'A';
        else if (a[i] == 'G') a[i] = 'C';
        else if (a[i] == 'C') a[i] = 'G';
    }

    reverse(a.begin(), a.end());
    return a;
}

BeAnalyzer::BeAnalyzer(const int R, const int filted_n, const int cleavage_length, const int pri_len, string grna_seq, string wt_seq) {
    int grna_pos, cleavage_position; // m_dir sequence direction

    if ((grna_pos = wt_seq.find(grna_seq)) != string::npos) {
        cleavage_position = grna_pos + cleavage_length;
        m_dir = 0;
    }
    else if ((grna_pos = wt_seq.find(rev_comp(grna_seq))) != string::npos) {
        cleavage_position = grna_pos + grna_seq.size() - cleavage_length;
        cleavage_position = wt_seq.size()-cleavage_position-1;
        wt_seq = rev_comp(wt_seq);
        grna_seq = rev_comp(grna_seq);
        m_dir = 1;
    }

    m_pri_for = wt_seq.substr(cleavage_position - R, pri_len);
    m_pri_back = wt_seq.substr(cleavage_position + R - pri_len, pri_len);

    m_wt_seq_sliced = wt_seq.substr(cleavage_position - R, 2 * R);

    m_cleavage_position = R;
    if (cleavage_position < R) m_cleavage_position -= R-cleavage_position;

    m_pattern = new sub_t[R * 2];
    m_cnt_position = new int[R * 2];

    m_cnt_pri = 0; m_cnt_filt = 0; m_cnt_all=0;
    m_pri_len = pri_len;
    m_R = R;
    m_grna_seq = grna_seq;
    return;
}

BeAnalyzer::~BeAnalyzer() {
    delete [] m_pattern;
    delete [] m_cnt_position;
}

int find_pri(string seq, string pri) {
    char ATGC[4] = { 'A','T','G','C' };
    int i, x, y, pos = -1;
    string mut_seq;

    for (y = 0; y < pri.size(); y++) {
        for (i = 0; i < 4; i++) {
            mut_seq = pri.substr(0, y) + ATGC[i] + pri.substr(y + 1, pri.size() - y - 1);
            if (seq.find(mut_seq) != string::npos) {
                pos = seq.find(mut_seq);
                break;
            }
        }
        if (pos >= 0) break;
    }
    return pos;
}

void BeAnalyzer::update_seq(string strseq) {
    int pos_for, pos_back, i = 0;
    string seq_sliced;

    if (m_dir == 1)
        strseq = rev_comp(strseq);

    m_cnt_all++;
    pos_for = find_pri(strseq, m_pri_for);
    pos_back = find_pri(strseq, m_pri_back);

    if (pos_for >= pos_back || pos_for < 0 || pos_back < 0) {
        return;
    }
    else {
        seq_sliced = strseq.substr(pos_for, pos_back - pos_for + m_pri_len);
        int comp = seq_sliced.find("N");
        if (comp > 0)
            return;
        m_cnt_pri++;
        m_cnt_filt++;
    }

    if (m_seq_count.find(seq_sliced) == m_seq_count.end()) {
        m_seq_count[seq_sliced] = 1;
    } else {
        m_seq_count[seq_sliced]++;
    }
}

void BeAnalyzer::read_fastq_files(const char* fastq_file1, const char* fastq_file2) {
    be_fastq_join(fastq_file1, fastq_file2);
}

void BeAnalyzer::read_fastq_files(const char* fastq_joined_file) {
    int linecnt;

    struct fq entry;
    meminit(entry);

    FILE* fp;
    fp = fopen(fastq_joined_file, "r");
    if (fp == NULL) {
        printf("can't find fastq file");
    }
    else {
        while (read_fq(fp, linecnt++, &entry)) {
            update_seq(string(entry.seq.s));
        }
    }
    fclose(fp);
}

void needle(string wt_seq, string sliced_seq, const int gapopen, const float gapextend, const int endgapopen, const float endgapextend, char *wt_str_char, char *seq_str_char, char *sym_str_char) {
    EM_ASM_({
        var rtn = needle(Pointer_stringify($0), Pointer_stringify($1), $2, $3, $4, $5);
        writeStringToMemory(rtn[0], $6);
        writeStringToMemory(rtn[1], $8);
        writeStringToMemory(rtn[2], $7);
    }, wt_seq.c_str(), sliced_seq.c_str(), gapopen, gapextend, endgapopen, endgapextend, wt_str_char, seq_str_char, sym_str_char);
}

void BeAnalyzer::run_alignment() {
    const int gapopen = 10;
    const float gapextend = 0.5;
    const int endgapopen = 10;
    const float endgapextend = 0.5;

    char *wt_str_char = new char[10000];
    char *seq_str_char = new char[10000];
    char *sym_str_char = new char[10000];

    m_sorted_list = sorted_filtered_map(m_seq_count, 1, &m_cnt_filt); //sort

    int prev_prog = 0;
    for (int i=0; i<m_sorted_list.size(); i++) {
        int new_prog = int((i/(float)m_sorted_list.size())*30);
        if (prev_prog < new_prog) {
            report_progress(i, m_sorted_list.size(), 30, 30, "Aligning reads...");
            prev_prog = new_prog;
        }

        needle(m_wt_seq_sliced, m_sorted_list[i], gapopen, gapextend, endgapopen, endgapextend, wt_str_char, seq_str_char, sym_str_char);

        m_emboss_wt_list.push_back(string(wt_str_char));
        m_emboss_seq_list.push_back(string(seq_str_char));
        m_emboss_sym_list.push_back(string(sym_str_char));
    }

    delete [] wt_str_char;
    delete [] seq_str_char;
    delete [] sym_str_char;
}

void BeAnalyzer::data_analyze() {
    m_cnt_insertion = 0;
    m_cnt_deletion = 0;
    m_cnt_wt = 0;
    m_cnt_sub = 0;
    m_cnt_others = 0;
    m_cnt_c_to_d = 0;


    int i = 0, j = 0, prev_prog = 0;

    bool has_c_to_d;
    for (vector<string>::iterator it=m_sorted_list.begin(); it!=m_sorted_list.end(); ++it) {
        int new_prog = int((i/(float)m_sorted_list.size())*30);
        if (prev_prog < new_prog) {
            report_progress(i, m_sorted_list.size(), 60, 30, "Analyzing results...");
            prev_prog = new_prog;
        }
        int this_cnt = m_seq_count[*it];
        if (m_emboss_wt_list[i].find('-') != string::npos && m_emboss_seq_list[i].find('-') != string::npos) {
            m_cnt_others += this_cnt;
            m_type.push_back(Other);
        }
        else if (m_emboss_wt_list[i].find('-') != string::npos) {
            m_cnt_insertion += this_cnt;
            m_type.push_back(Ins);
        }
        else if (m_emboss_seq_list[i].find('-') != string::npos) {
            m_cnt_deletion += this_cnt;
            m_type.push_back(Del);
        }
        else {
            if (m_emboss_sym_list[i].find('.') != string::npos) {
                m_cnt_sub += this_cnt;
                m_type.push_back(Sub);
                has_c_to_d = false;
                for (j = 0; j < m_R * 2; j++) {
                    if (m_emboss_sym_list[i].substr(j, 1) == ".") {
                        if (m_emboss_wt_list[i].substr(j, 1) == "C") {
                            if (!has_c_to_d) {
                                m_cnt_c_to_d += this_cnt;
                                has_c_to_d = true;
                            }
                        }
                        m_cnt_position[j] += this_cnt;
                        m_pattern[j][m_emboss_seq_list[i][j]] += this_cnt;
                    }
                }
            }
            else {
                m_type.push_back(WT);
                m_cnt_wt += this_cnt;
            }
        }
        i++;
    }
    return;
}

void BeAnalyzer::write_align() {
    int i;
    FILE* fp = stdout;
    vector<string> types = {"WT", "Ins", "Del", "Sub", "Other"};

    for (i = 0; i < m_sorted_list.size(); i++) {
        fprintf(fp, "%s\t%s\t%d\n%s\n%s\n\n", m_emboss_wt_list[i].c_str(), types[m_type[i]].c_str(), m_seq_count[m_sorted_list[i]], m_emboss_sym_list[i].c_str(), m_emboss_seq_list[i].c_str());
    }
    fclose(fp);
    return;
}

void BeAnalyzer::write_count() {

    //FILE* fp = fopen("count.txt", "w");
    FILE *fp = stdout;
    fprintf(fp, "All\t%d\nPrimer\t%d\nfilted\t%d\nWild type\t%d\ninsertion\t%d\ndeletion\t%d\nsubstitution\t%d\nC to D\t%d\n", m_cnt_all, m_cnt_pri, m_cnt_filt, m_cnt_wt, m_cnt_insertion, m_cnt_deletion, m_cnt_sub, m_cnt_c_to_d);
    fclose(fp);
    return;
}

void BeAnalyzer::write_substitution() {

    int i = 0, c = 0;

    //FILE* fp = fopen("/home/parkj/experiments/Be_Analyzer/Total_Substitution.txt", "w");
    FILE* fp = stdout;
    fprintf(fp, "Site\tWt_seq\tSubstitution_count\tA\tT\tG\tC\n");
    c = 0;
    int loop_start = -m_cleavage_position;
    int loop_end = m_wt_seq_sliced.size() - m_cleavage_position + 1;
    for (i = loop_start; i < loop_end; i++) {
        if (i == 0) {
            continue;
        }
        else {
            c += 1;
            fprintf(fp, "%d\t%s\t%d\t%d\t%d\t%d\t%d\n", i, m_wt_seq_sliced.substr(c - 1, 1).c_str(), m_cnt_position[c - 1], m_pattern[c - 1]['A'], m_pattern[c - 1]['T'], m_pattern[c - 1]['G'], m_pattern[c - 1]['C']);
        }
    }
    fclose(fp);
    return;
}
