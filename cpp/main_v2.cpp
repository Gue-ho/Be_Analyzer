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
#include "fastq-lib.h"

using namespace std;

enum BeType {
    WT = 0,
    Ins = 1,
    Del = 2,
    Sub = 3,
    Other = 4
};

const int R = 60, filted_n = 1, cleavage_length = 17, pri_len = 10; // 17bp + cleavage + 3bp + NGG for cas9
string grna_seq = "GAGTCCGAGCAGAAGAAGAA"; // 20bp for cas9
string wt_seq = "GGGCCTCCTGAGTTTCTCATCTGTGCCCCTCCCTCCCTGGCCCAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAACCACAAACCCACGAGGGCAGAGTGCTGCTTGCTGCTG"; // wt_seq for cas9
const char* fastq_joined_file = "/home/parkj/experiments/Be_Analyzer/data/7.fastqjoin"; // index.fastqjoin
string control_file = "/home/parkj/experiments/Be_Analyzer/data/15.fastqjoin";

typedef map<string, int> seq_count;

struct sub_t {
    int cnt_A;
    int cnt_T;
    int cnt_G;
    int cnt_C;
    int& operator[](char sub) {
        if (sub == 'A') return cnt_A;
        if (sub == 'T') return cnt_T;
        if (sub == 'G') return cnt_G;
        if (sub == 'C') return cnt_C;
    }
    const int& operator[](char sub) const {
        if (sub == 'A') return cnt_A;
        if (sub == 'T') return cnt_T;
        if (sub == 'G') return cnt_G;
        if (sub == 'C') return cnt_C;
    }
};


struct less_second {
    typedef pair<string, int> type;
    bool operator () (type const& a, type const &b) const {
        return a.second < b.second;
    }
};
vector<string> sorted_map(seq_count m_out_seq) {
	vector<pair<string, int> > cp(m_out_seq.begin(), m_out_seq.end());
	sort(cp.begin(), cp.end(), less_second());

    vector<string> rtn;
    for (int i=0; i<cp.size(); i++)
        rtn.push_back(cp[i].first);

	return rtn;
}

class BeAnalyzer {
public:
    ~BeAnalyzer();
	BeAnalyzer();
	void process();
	void data_analyze();
	void write_align();
	void write_count();
	void write_substitution();
	void write_c_to_d();

	int m_dir, m_cnt_all, m_cnt_pri, m_cnt_filt;
	string m_wt_seq_sliced, m_pri_for, m_pri_back;
	int m_cnt_insertion, m_cnt_deletion, m_cnt_wt, m_cnt_others;
	string m_position, m_c_to_d_position;

    vector<string> m_emboss_wt_list, m_emboss_seq_list, m_emboss_sym_list, m_sorted_list;

	seq_count m_seq_count;
    sub_t* m_pattern;
    sub_t* m_c_to_d_sub;

    int *m_cnt_sub;
    int *m_cnt_c_to_d;

    vector<int> m_type;
};

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

int find_string(string grans_seq) {
	if (wt_seq.find(grans_seq) != string::npos) return wt_seq.find(grans_seq);
	else if (wt_seq.find(rev_comp(grans_seq)) != string::npos) return wt_seq.find(rev_comp(grans_seq));
}

BeAnalyzer::BeAnalyzer() {
	int grna_pos, cleavage_position; // m_dir sequence direction

	if ((grna_pos = wt_seq.find(grna_seq)) != string::npos) {
		cleavage_position = grna_pos - cleavage_length;
	}
	else if ((grna_pos = wt_seq.find(rev_comp(grna_seq))) != string::npos) {
		cleavage_position = grna_pos - grna_seq.size() + cleavage_length;
		m_dir = 1;
	}

	m_pri_for = wt_seq.substr(cleavage_position - R, pri_len);
	m_pri_back = wt_seq.substr(cleavage_position + R - pri_len, pri_len);

	m_wt_seq_sliced = wt_seq.substr(cleavage_position - R, 2 * R);

    m_pattern = new sub_t[R * 2];
    m_c_to_d_sub = new sub_t[R * 2];

    m_cnt_sub = new int[R * 2];
    m_cnt_c_to_d = new int[R * 2];

	return;
}


BeAnalyzer::~BeAnalyzer() {
    delete [] m_pattern;
    delete [] m_c_to_d_sub;

    delete [] m_cnt_sub;
    delete [] m_cnt_c_to_d;
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

void BeAnalyzer::process() {

	int line = 0, pos_for, pos_back, i = 0;
    m_cnt_pri = 0; m_cnt_filt = 0; m_cnt_all=0;
	string strseq, seq_sliced;

    struct fq entry;
    meminit(entry);

	FILE* fp;
	fp = fopen(fastq_joined_file, "r");
	if (fp == NULL) {
		printf("can't find fastq file");
	}
	else {
        while (read_fq(fp, line, &entry)) {
			strseq = string(entry.seq.s);
			m_cnt_all++;
			pos_for = find_pri(strseq, m_pri_for);
			pos_back = find_pri(strseq, m_pri_back);

			if (pos_for >= pos_back || pos_for < 0 || pos_back < 0) {
				continue;
			}
			else {
				m_cnt_pri++;
				m_cnt_filt++;
				seq_sliced = strseq.substr(pos_for, pos_back - pos_for + pri_len);
			}

			if (m_seq_count.find(seq_sliced) == m_seq_count.end())
				m_seq_count[seq_sliced] = 1;
			else
        		m_seq_count[seq_sliced]++;

            line++;
		}
	}
	fclose(fp);

	for (seq_count::iterator iter = m_seq_count.begin(); iter != m_seq_count.end(); ++iter) {
		if (iter->second == 1) {
			m_cnt_filt--;
			m_seq_count.erase(iter);
			continue;
		}
	}
    m_sorted_list = sorted_map(m_seq_count); //sort

	fp = fopen("/home/parkj/experiments/Be_Analyzer/emboss_wt.txt", "w");
	fprintf(fp, ">\n%s\n", m_wt_seq_sliced.c_str());
	fclose(fp);

	fp = fopen("/home/parkj/experiments/Be_Analyzer/emboss_before.txt", "w");
	for (i=0; i<m_sorted_list.size(); i++) {
		fprintf(fp, ">\n%s\n", m_sorted_list[i].c_str());
	}
	fclose(fp);

	system("needle /home/parkj/experiments/Be_Analyzer/emboss_wt.txt /home/parkj/experiments/Be_Analyzer/emboss_before.txt -gapopen 10 -gapextend 0.5 -outfile /home/parkj/experiments/Be_Analyzer/emboss_after.txt");
	return;
}

void BeAnalyzer::data_analyze() {
	m_cnt_insertion = 0;
	m_cnt_deletion = 0;
	m_cnt_wt = 0;
	m_cnt_sub = 0;
	m_cnt_others = 0;
    m_cnt_c_to_d = 0;

	int c = 0, w = 0, i = 0;
	string emboss_wt = "", emboss_seq = "", emboss_sym = "", line;
	char line_char[500];

	FILE* fp = fopen("emboss_result.txt", "r");
	while (fgets(line_char, 500, fp)) {
		line = (string) line_char;
		if (line.find("EMBOSS_001") == string::npos) continue;
		else {
			c = 1;
			emboss_wt += line.substr(21, 50);
			w += 1;
		}
		if (w != 0) w++;
		if (w == 3) emboss_sym += line.substr(21, 50);
		if (w == 4) {
			w = 0;
			emboss_seq += line.substr(21, 50);
		}
		if (c == 1 && line.find("====") != string::npos) {
			m_emboss_wt_list.push_back(emboss_wt);
			m_emboss_seq_list.push_back(emboss_seq);
			m_emboss_sym_list.push_back(emboss_sym);
			emboss_wt = ""; emboss_seq = ""; emboss_sym = "";
		}
	}
	fclose(fp);

	for (seq_count::iterator it=m_seq_count.begin(); it!=m_seq_count.end(); ++it) {
		if (m_emboss_wt_list[i].find('-') != string::npos) {
			m_cnt_others += it->second;
			m_type[i] = Other;
		}
		else if (m_emboss_wt_list[i].find('-') != string::npos) {
			m_cnt_insertion += it->second;
			m_type[i] = Ins;
		}
		else if (m_emboss_seq_list[i].find('-') != string::npos) {
			m_cnt_deletion += it->second;
			m_type[i] = Del;
		}
		else {
			if (m_emboss_sym_list[i].find('.') != string::npos) {
				m_cnt_sub += it->second;
				m_type[i] = Sub;
				for (i = 0; i < R * 2; i++) {
					if (m_emboss_sym_list[i].substr(i, 1) == ".") {
						if (m_emboss_wt_list[i].substr(i, 1) == "C") {
                            m_cnt_c_to_d += it->second;
							m_c_to_d_position[i] += it->second;
							m_c_to_d_sub[i][m_emboss_wt_list[i][0]] += it->second;
						}
						m_position[i] += it->second;
						m_pattern[i][m_emboss_wt_list[i][0]] += it->second;
					}
				}
			}
			else {
                m_type[i] = WT;
				m_cnt_wt += it->second;
			}
		}
	}
	return;
}

void BeAnalyzer::write_align() {
	int i;

	FILE* fp = fopen("align_result.txt", "w");
	for (i = 0; i < m_sorted_list.size(); i++) {
		fprintf(fp, "%s\t%s\t%d\n%s\n%s\n\n", m_emboss_wt_list[i].c_str(), m_type[i], m_seq_count[m_sorted_list[i]], m_emboss_sym_list[i].c_str(), m_emboss_seq_list[i].c_str());
	}
	fclose(fp);
	return;
}

void BeAnalyzer::write_count() {

	FILE* fp = fopen("count.txt", "w");
	fprintf(fp, "All\t%d\nPrimer\t%d\nfilted\t%d\nWild type\t%d\ninsertion\t%d\ndeletion\t%d\nsubstitution\t%d\nC to D\t%d\n", m_cnt_all, m_cnt_pri, m_cnt_filt, m_cnt_wt, m_cnt_insertion, m_cnt_deletion, m_cnt_sub, m_cnt_c_to_d);
	fclose(fp);
	return;
}

void BeAnalyzer::write_substitution() {

	int i = 0, c = 0;

	FILE* fp = fopen("Total_Substitution.txt", "w");
	fprintf(fp, "Site\tWt_seq\tSubstitution_count\tA\tT\tG\tC\n");
	c = 0;
	for (i = (m_wt_seq_sliced.find(grna_seq) + 17) * (-1); i < (m_wt_seq_sliced.find(grna_seq) + 17) * (-1) + m_wt_seq_sliced.size() + 1; i++) {
		if (i == 0) {
			continue;
		}
		else {
			c += 1;
			fprintf(fp, "%d\t%s\t%d\t%d\t%d\t%d\t%d\n", i, m_wt_seq_sliced.substr(c - 1, 1).c_str(), m_position[c - 1], m_pattern[c - 1]['A'], m_pattern[c - 1]['T'], m_pattern[c - 1]['G'], m_pattern[c - 1]['C']);
		}
	}
	fclose(fp);
	return;
}

void BeAnalyzer::write_c_to_d() {
	int i = 0, c = 0;

	FILE* fp = fopen("C_to_D_Substituion.txt", "w");
	fprintf(fp, "Site\tWt_seq\tSubstitution_count\tA\tT\tG\tC\n");
	c = 0;
	for (i = (m_wt_seq_sliced.find(grna_seq) + 17) * (-1); i < (m_wt_seq_sliced.find(grna_seq) + 17) * (-1) + m_wt_seq_sliced.size() + 1; i++) {
		if (i == 0) {
			continue;
		}
		else {
			c += 1;
			fprintf(fp, "%d\t%s\t%d\t%d\t%d\t%d\t%d\n", i, m_wt_seq_sliced.substr(c - 1, 1).c_str(), m_c_to_d_position[c - 1], m_c_to_d_sub[c - 1]['A'], m_c_to_d_sub[c - 1]['T'], m_c_to_d_sub[c - 1]['G'], m_c_to_d_sub[c - 1]['C']);
		}
	}
	fclose(fp);
	return;
}

int main (void) {
    BeAnalyzer ba = BeAnalyzer();
    ba.process();
    puts("Hello World!");
    return 0;
}