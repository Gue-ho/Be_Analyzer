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

int R = 60, filted_n = 1, cleavage_length = 17; // 17bp + cleavage + 3bp + NGG for cas9
std::string grna_seq = "ATAGATAGATGTAGAGT"; // 20bp for cas9
std::string wt_seq = ""; // wt_seq for cas9
std::string fastq_joined_file = "7.fastqjoin"; // index.fastqjoin
std::string control_file = "15.fastqjoin";

typedef std::map<std::string, int> seq_count;
typedef std::map<int, std::vector<std::string>> pattern;

class BeAnalyzer {
private:
	BeAnalyzer();
	void process(std::string, int);
	void data_analyze(std::string, std::string);
	void write_align(std::string, std::string, std::string, std::string, std::string, std::string);
	void write_count(int, int, int, int, int, int, int, int);
	void write_substitution(std::string, std::string, std::string, pattern);
	void write_c_to_d(std::string, std::string, std::string, pattern);

	int m_dir, cnt_all, cnt_pri, cnt_filt;
	std::string m_wt_seq_sliced, m_pri_for, m_pri_back, m_sorted_list, m_seq_count;
	int cnt_all, cnt_pri, cnt_filt, cnt_insertion, cnt_deletion, cnt_wt, cnt_sub, cnt_c_to_d, cnt_others;
	std::string m_position, m_c_to_d_position, m_emboss_wt_list, m_emboss_seq_list, m_emboss_sym_list, m_ori_seq_list, m_type;

	seq_count m_seq_count;
	pattern m_pattern, m_c_to_d_sub;
};

std::string rev_comp(std::string &a) {
	//std::string a must be uppercase letters
	int i;
	std::string s;

	for (i = 0; i<a.size(); i++) {
		if (a[i] == 'A') a[i] = 'T';
		else if (a[i] == 'T') a[i] = 'A';
		else if (a[i] == 'G') a[i] = 'C';
		else if (a[i] == 'C') a[i] = 'G';
	}

	reverse(a.begin(), a.end());
	return a;
}

int find_string(std::string grans_seq) {
	if (wt_seq.find(grans_seq) != std::string::npos) return wt_seq.find(grans_seq);
	else if (wt_seq.find(rev_comp(grans_seq)) != std::string::npos) return wt_seq.find(rev_comp(grans_seq));
}

BeAnalyzer::BeAnalyzer() {
	int grna_pos, cleavage_position, m_dir = 0; // m_dir sequence direction
	std::string m_pri_for, m_pri_back, m_wt_seq_sliced;

	if (wt_seq.find(grna_seq) != std::string::npos) {
		grna_pos = wt_seq.find(grna_seq);
		cleavage_position = grna_pos - cleavage_length;
	}
	else if (wt_seq.find(grna_seq) != std::string::npos) {
		grna_pos = wt_seq.find(rev_comp(grna_seq));
		cleavage_position = grna_pos - grna_seq.size() + cleavage_length;
		wt_seq = rev_comp(wt_seq);
		m_dir = 1;
	}

	m_pri_for = wt_seq.substr(cleavage_position - R, 10);
	m_pri_back = wt_seq.substr(cleavage_position + R - 10, 10);

	m_wt_seq_sliced = wt_seq.substr(cleavage_position - R, 2 * R);

	return;
}

int find_pri(std::string seq, std::string pri) {
	char ATGC[4] = { 'A','T','G','C' };
	int i, x, y, pos = 0;
	std::string mut_seq;

	for (y = 0; y < pri.size(); y++) {
		for (i = 0; i < 4; i++) {
			mut_seq = pri.substr(0, y) + ATGC[i] + pri.substr(y + 1, 10 - y - 1);
			if (seq.find(mut_seq) != std::string::npos) {
				pos = seq.find(mut_seq);
				break;
			}
		}
	}
	return pos;
}

void BeAnalyzer::process(std::string m_wt_seq_sliced, int m_dir) {

	int line = 0, cnt_all = 0, cnt_pri = 0, cnt_filt = 0, pos_for, pos_back, i = 0;
	typedef std::map <std::string, int> seq_count;
	typedef std::map <int, std::string> seq_list;

	char seq[500];
	std::string strseq, seq_sliced;

	seq_count m_seq_count;
	seq_list seq_list_map;

	FILE* fp;
	fopen_s(&fp, fastq_joined_file, "r");
	if (fp == NULL) {
		printf("can't find fastq file");
	}
	else {
		while (fgets(seq, 500, fp)) {
			line++;
			strseq = (std::string) seq;
			if (line == 1) {
				cnt_all++;
				line = 0;
				pos_for = find_pri(seq, m_pri_for);
				pos_back = find_pri(seq, m_pri_back);

				if (pos_for >= pos_back) {
					continue;
				}
				else {
					cnt_pri++;
					cnt_filt++;
					seq_sliced = strseq.substr(pos_for, pos_back * 2 - pos_for);
				}

				if (m_seq_count.find(seq_sliced) == m_seq_count.end()) {
					seq_list_map[i] = seq_sliced;
					i++;
					m_seq_count[seq_sliced] = 1;
				}
				else {
					m_seq_count[seq_sliced] ++;
				}
			}
		}
	}
	fclose(fp);

	int ii = 0;

	for (i = 0; i < m_seq_count.size(); i++) { //filter count 1
		if (m_seq_count[seq_list_map[i]] == 1) {
			cnt_filt += -1;
			m_seq_count.erase(seq_list_map[i]);
		}

		seq_list_map.clear();

		std::string m_sorted_list[m_seq_count.size()]; //sort
		int i = 0;

		typedef std::function < bool(std::pair<std::string, int>, std::pair <std::string, int>)> Comparator;
		Comparator compFunctor = [](std::pair<std::string, int> elem1, std::pair <std::string, int> elem2) {
			return elem1.second > elem2.second;
		};
		std::set <std::pair <std::string, int>, Comparator> setOfWords(m_seq_count.begin(), m_seq_count.end(), compFunctor);
		for (std::pair<std::string, int> element : setOfWords) {
			m_sorted_list[i] = element.first;
		}

		FILE* fp;
		fopen_s(&fp, "emboss_before.txt", "w");
		for (i = 0; i < m_sorted_list.size(); i++) {
			fprintf(fp, ">\n%s\n", m_sorted_list[i]);
		}
		fclose(fp);

		FILE* fp;
		fopen_s(&fp, "emboss_wt.txt", "w");
		fprintf(fp, ">\n%s\n", m_wt_seq_sliced);
		fclose(fp);

		system("needle emboss_wt.txt emboss_before.txt -gapopen 10 -gapextend 0.5 -outfile emboss_after.txt");

		return;
	}
}

void BeAnalyzer::data_analyze(std::string m_seq_count, std::string m_sorted_list) {
	int cnt_insertion = 0, cnt_deletion = 0, cnt_wt = 0, cnt_sub = 0, cnt_others = 0, cnt_c_to_d_sub = 0, c = 0, w = 0, i = 0;
	std::string emboss_wt = "", emboss_seq = "", emboss_sym = "", ori_seq = "", m_emboss_wt_list[m_sorted_list.size()], m_emboss_seq_list[m_sorted_list.size()], m_emboss_sym_list[m_sorted_list.size()], m_ori_seq_list[m_sorted_list.size()], line;
	char line_char[500];

	FILE* fp;
	fopen_s(&fp, "emboss_result.txt", "r");
	while (fgets(line_char, 500, fp)) {
		line = (std::string) line_char;
		if (line.find("EMBOSS_001") == std::string::npos) continue;
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
		if (c == 1 && line.find("====") != std::string::npos) {
			ori_seq = emboss_seq;
			while (true) {
				if (ori_seq.find('-') == std::string::npos) break;
				else ori_seq.erase(ori_seq.find('-', 1));
			}
			m_emboss_wt_list[i] = emboss_wt;
			m_emboss_seq_list[i] = emboss_seq;
			m_emboss_sym_list[i] = emboss_sym;
			m_ori_seq_list[i] = ori_seq;
			emboss_wt = "", emboss_seq = "", emboss_sym = "", ori_seq = "";
		}
	}
	fclose(fp);

	std::string m_pattern[4][R * 2], m_c_to_d_sub[4][R * 2], m_position[R * 2], m_c_to_d_position[R * 2], m_type[m_sorted_list.size()], ATGC = "ATGC";

	for (i = 0; i < m_ori_seq_list.sieq(), i++;) {
		if (m_emboss_wt_list[i].find('-') != std::string::npos) {
			cnt_others += m_seq_count[m_ori_seq_list[i]];
			m_type[i] = "others";
		}
		else if (m_emboss_wt_list[i].find('-') != std::string::npos) {
			cnt_insertion += m_seq_count[m_ori_seq_list[i]];
			m_type[i] = "insertion";
		}
		else if (m_emboss_seq_list[i].find('-') != std::string::npos) {
			cnt_deletion += m_seq_count[m_ori_seq_list[i]];
			m_type[i] = "deletion";
		}
		else {
			if (m_emboss_sym_list[i].find('.') != std::string::npos) {
				cnt_sub += m_seq_count[m_ori_seq_list[i]];
				m_type[i] = "Substitution";
				for (i = 0; i < R * 2; i++) {
					if (m_emboss_sym_list[i].substr(i, 1) == ".") {
						if (m_emboss_wt_list[i].substr(i, 1) == "C") {
							cnt_c_to_d_sub += m_seq_count[m_ori_seq_list[i]];
							m_c_to_d_position[i] += m_seq_count[m_ori_seq_list[i]];
							m_c_to_d_sub[ATGC.find(m_emboss_wt_list[i])][i] += m_seq_count[m_ori_seq_list[i]];
						}
						m_position[i] += m_seq_count[m_ori_seq_list[i]];
						m_pattern[ATGC.find(m_emboss_wt_list[i])][i] += m_seq_count[m_ori_seq_list[i]];
					}
				}
			}
			else {
				cnt_wt += m_seq_count[m_ori_seq_list[i]];
			}
		}
	}
	return;
}

void write_align(std::string m_emboss_wt_list, std::string m_emboss_seq_list, std::string m_emboss_sym_list, std::string m_seq_count, std::string m_sorted_list, std::string m_type) {

	int i;

	FILE* fp;
	fopen_s(&fp, "align_result.txt", "w");
	for (i = 0; i < m_sorted_list.size(); i++) {
		fprintf(fp, "%s\t%s\t%d\n%s\n%s\n\n", m_emboss_wt_list[i], m_type[i], m_seq_count[m_sorted_list[i]], m_emboss_sym_list[i], m_emboss_seq_list[i]);
	}
	fclose(fp);
	return;
}

void write_count(int cnt_all, int cnt_pri, int cnt_filt, int cnt_wt, int cnt_insertion, int cnt_deletion, int cnt_sub, int cnt_c_to_d_sub) {

	FILE* fp;
	fopen_s(&fp, "count.txt", "w");
	fprintf(fp, "All\t%d\nPrimer\t%d\nfilted\t%d\nWild type\t%d\ninsertion\t%d\ndeletion\t%d\nsubstitution\t%d\nC to D\t%d\n", cnt_all, cnt_pri, cnt_filt, cnt_wt, cnt_insertion, cnt_deletion, cnt_sub, cnt_c_to_d_sub);
	fclose(fp);
	return;
}

void write_substitution(std::string m_wt_seq_sliced, std::string grna_seq, std::string m_position, pattern m_pattern) {

	int i = 0, c = 0;

	FILE* fp;
	fopen_s(&fp, "Total_Substitution.txt", "w");
	fprintf(fp, "Site\tWt_seq\tSubstitution_count\tA\tT\tG\tC\n");
	c = 0;
	for (i = (m_wt_seq_sliced.find(grna_seq) + 17) * (-1); i < (m_wt_seq_sliced.find(grna_seq) + 17) * (-1) + m_wt_seq_sliced.size() + 1; i++) {
		if (i == 0) {
			continue;
		}
		else {
			c += 1;
			fprintf(fp, "%d\t%s\t%d\t%d\t%d\t%d\t%d\n", i, m_wt_seq_sliced.substr(c - 1, 1), m_position[c - 1], m_pattern['A'][c - 1], m_pattern['T'][c - 1], m_pattern['G'][c - 1], m_pattern['C'][c - 1]);
		}
	}
	fclose(fp);
	return;
}

void write_c_to_d(std::string m_wt_seq_sliced, std::string grna_seq, std::string m_c_to_d_position, pattern m_c_to_d_sub) {
	int i = 0, c = 0;

	FILE* fp;
	fopen_s(&fp, "C_to_D_Substituion.txt", "w");
	fprintf(fp, "Site\tWt_seq\tSubstitution_count\tA\tT\tG\tC\n");
	c = 0;
	for (i = (m_wt_seq_sliced.find(grna_seq) + 17) * (-1); i < (m_wt_seq_sliced.find(grna_seq) + 17) * (-1) + m_wt_seq_sliced.size() + 1; i++) {
		if (i == 0) {
			continue;
		}
		else {
			c += 1;
			fprintf(fp, "%d\t%s\t%d\t%d\t%d\t%d\t%d\n", i, m_wt_seq_sliced.substr(c - 1, 1), m_c_to_d_position[c - 1], m_c_to_d_sub['A'][c - 1], m_c_to_d_sub['T'][c - 1], m_c_to_d_sub['G'][c - 1], m_c_to_d_sub['C'][c - 1]);
		}
	}
	fclose(fp);
	return;
}

