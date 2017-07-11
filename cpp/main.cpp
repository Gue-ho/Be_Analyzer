#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <math.h>
#include <ctype.h>

using namespace std;

typedef map<String, String> emboss;
typedef map<String, int> Sequence;
class BeAnalyzer {
    private:
        // String emobossfile='emboss_after.txt' //emboss result file name
        BeAnalyzer();
        void seq_standard(String, String, int, int);
        void indicator_table(String, int);
        void read_seq(String, String, String, int);
        void read_EMBOSS(String, Sequence, String);
        void sort(Sequence);
        void write_emboss_dict(String, Sequence);
        void write_emboss_result(String, Sequence);
        void write_total_substitution(String, String, String, String);
        void write_C_to_D_substitution(String, String, String);

        int m_all_count, m_primer_count, m_filted_count;
        int m_count_list[6];

        Sequence m_out_seq;

        String m_s_seq, m_seq_range;
        String m_for_table, m_back-table;

        vector<String> m_substitution_position;
        vector<String> m_sorted_list;

        map<char, vector<String> > m_substitution_pattern;
        map<char, vector<String> > m_C_to_D;

        emboss m_emboss_dict;
}

void rev_comp(String &a) {
    // String a must be uppercase letters
    int i;
    for (i=0; i<a.size(); i++) {
        if (a[i] == 'A') a[i] = 'T';
        else if (a[i] == 'T') a[i] = 'A';
        else if (a[i] == 'G') a[i] = 'C';
        else if (a[i] == 'C') a[i] = 'G';
    }

    reverse(a.begin(), a.end());
}

void BeAnalyzer::BeAnalyzer() {
    // do nothing
    return;
}

void BeAnalyzer::seq_standard(String wt_seq, String target_seq, int filt_r, int end_range) {
    int i, start_pos, end_pos;

    for (i=22; i<wt_seq.size()-23; i++) {  //Todo:SpCas9 target seq
        if (wt_seq.substr(i, 23).compare(target_seq) == 0) {
            start_pos = i+ 17 - end_range;
            end_pos = i + 17 + end_range; //Todo:SpCas9 celvage
            
            if (start_pos < 0) start_pos = 0;
            if (end_pos > wt_seq.size()) end_pos = wt_seq.size();

            m_s_seq = wt_seq.substr(i+17-filt_r, filt_r*2);
            m_seq_range = wt_seq.substr(start_pos, end_pos - start_pos + 1);

            break;
        } else if (wt_seq.substr(i, 23).compare(rev_comp(target_seq))) {
            start_pos = 1 + 6 - end_range;
            end_pos = i+6 + end_range;

            if (start_pos < 0) start_pos = 0;
            if (end_pos > wt_seq.size()) end_pos = wt_seq.size();

            m_s_seq = wt_seq.substr(i+6-filt_r, filt_r*2);
            m_seq_range = wt_seq.substr(start_pos, end_pos - start_pos + 1);
            break
        } else {
            start_pos = 0;
            end_pos = wt_seq.size();

            m_s_seq = "None";
            m_seq_range = wt_seq;
        }
    }
    return;
}

void BeAnalyzer::indicator_table(int indicator_range) {
    String pri_for=m_seq_range.substr(0,indicator_range);
    String pri_back=m_seq_range.substr(indicator_range.size()-indicator_range,indicator_range);
    length_range=m_seq_range.size();
    String m_for_table[pow(4,indicator_range)],m_back_table[pow(4,indicator_range)];
    int i, x;
    String nucleotide[]={'A','T','G','C'};
    for (i,i<indicator_range.size();i++){
        for (x;x<4,x++){
            m_for_table[i]=pri_for.substr(0,i)+nucleotide[x]+pri_for.substr(i,indicatro_range);
            m_back_table[i]=pri_back.substr(0,i)+nucleotide[x]+pri_back.substr(i,indicator_range);
        }
    }
    return;
}

const String seq_check(String out_seq, String for_table, String back_table) {
    int for_check=0, back_check=0, i, start_seq, end_seq;
    String N='N', acu_seq;
    if (out_seq.find(N)!=string::npos){
        return '0';
    }
    
    for (i, i<for_table.size(),i++){
        if (out_seq.find(for_table[i]!=-1) {
            for_check=1;
            start_seq=out_seq.find(for_table[i]);
        }
        if (out_seq.find(back_table[i]!=-1) {
            back_check=1;
            back_seq=out_seq.find(back_table[i]);
        }

    if (for_check==1 && back_check==1) {
        acu_seq=out_seq.substr(start_seq, start_seq - end_seq + 1);
        return acu_seq;
    } else {
        return '0';
}

/*
class Sequence {
    public:
        String seq; //end_range pudlic?
        int score;
};
*/

typedef map<String, int> Sequence;
typedef map<int, String> Filt;
    
void BeAnalyzer::read_seq(String file_name, m_for_table, m_back_table, filt_n) {
    int m_all_count=0, m_primer_count=0, m_filted_count=0, line=0, read=0, i, h, fread=0;
    String* seq;
    const String acu_seq;

    Sequence m_out_seq;
    Filt filted_seq;

    FILE* fp;
    fopen_s(&fp, file_name,'r');
    if (fp==NULL); {
        printf('can't find fastq file');
    } else {
        while (fgets(seq,500,fp)) {
            line+=1;
            h=0;
            if (line==1) {
                line=0;
                m_all_count+=1;
                acu_seq=seq_check(sesq, m_for_table, m_back_table);
                if (acu_seq==0) {
                    continue;
                }
                if (out_seq.find(acu_seq) == out_seq.end()){
                    read+=1;
                    m_out_seq[acu_seq] = 1;
                    filted_seq[fread] = acu_seq;
                    fread+=1;
                } else {
                    m_out_seq[acu_seq]++;   
                }
            }
        }            
    }
    fclose(fp);
    for (i=0, i<fread, i++) {
        if (m_out_seq[filted_seq[i]]<=filted_n);
            m_out_seq.erase(filted_seq[i]);
            filted_count+=1
            continue;
        }
    return;
    }
}

/*
void emboss_needle () {
}
*/

int bar_count(String a){
    int c=0;
    while (true) {
        if (a.find('-')==str::string::npos) {
            break;
        } else {
            c+=1;
            a.earse(a.find('-'),1);
        }
    return c;
}

String str_lower(String a) {
    if (a=='A') a='a';
    else if (a=='T') a='t';
    else if (a=='G') a='g';
    else if (a=='C') a='c';
    return a;
}
void BeAnalyzer::read_EMBOSS(String file_name, Sequence m_out_seq, String m_seq_range) {
    
    int i, w=0, c=0;
    int m_count_list[6]={0}; //[WT,Substitution,CtoD,INs,Del,Others]
    String m_substitution_position[m_seq_range.size()]={0}, m_substitution_pattern[4][m_seq_range.size()]={0}, m_C_to_D[4][m_seq_range.size()]={0}, emboss_wt='', emboss_seq='', line, ori_seq;
    typedef map <String, String> emboss;
    emboss m_emboss_dict;
    
    FILE* fp;
    fopen(&fp, file_name, 'r');
    if (fp==NULL) {
        printf('can't find emboss result');
    } else {
        while (fgets(line, 500, fp)) {
            if (line.find('EMBOSS_001')==std::string::npos) {
                continue;
            } else {
                c=1;
                emboss_wt+=line.substr(line.find('EMBOSS_001')+21,50);
                w+=1;
            }
            if (w!=0) {
                w+=1;
            }
            if (w==4) {
                w=0;
                emboss_seq+=line.substr(line.find('EMBOSS'))+21,50);
            }
            if (c==1 && line.find('====')!=std::string::npos) {
                ori_seq=emboss_seq;
                while (true) {
                    if (ori_seq.find('-')==std::string::npos) {
                        break;
                    } else {
                        ori_seq.erase(ori_seq.find('-'),1);
                    }
                }
            
                String sub_wt='', sub_seq='', sub_sym=''; //Substitution
                
                for (i=0, i<emboss_wt.size(), i++) {
                    if (emboss_wt.substr(i,1)!='-' && emboss_seq.substr(i,1)!='-') {
                        if (emboss_wt.substr(i,1)!=emboss_seq.substr(i,1)) {
                            sub_wt+=emboss_wt.substr(i);
                            sub_seq+=str_lower(emboss_seq.substr(i));
                            sub_sym+='.';
                        } else {
                            sub_wt+=emboss_wt.substr(i,1);
                            sub_seq+=emboss_seq.substr(i,1);
                            sub_sym+='|';
                        }
                    } else {
                        sub_wt+=emboss_wt.substr(i,1);
                        sub_seq+=emboss_seq.substr(i,1);
                        sub_sym+=' ';
                    }
                }
                // Indel count
                String indel, atgc='atgc';
                int insbar=bar_count(emboss_wt), delbar=bar_count(emboss_seq);
                if (emboss_wt.find('-')!=str::string::npos && emboss_seq.find('-')!-str::string::npos) {
                    indel='Others';
                    m_count_list[5]+=m_out_seq[ori_seq];
                } else if (insbar>delbar) {
                    indel='Insertion';
                    m_count_list[3]+=m_out_seq[ori_seq];
                } else if (delbar>insbar) {
                    indel='Deletion';
                    m_count_list[4]+=m_out_seq[ori_seq];
                } else {
                    if (sub_sym.find('.')!=str::string::npos) {
                        indel='substitution';
                        m_count_list[1]+=m_out_seq[ori_seq];
                        for (i=0,i<sub_seq.size(),i++) {
                            if (atgc.find(sub_seq.substr(i,1))!=str::string::npos) {
                                substitution_position[i]+=m_out_seq[ori_seq]
                                substitution_pattern[atgc.find(sub_seq.substr(i,1))][i]+=m_out_seq[ori_seq];
                                if (sub_wt.strsub(i,1)=='C') {
                                    m_count_list[2]+=m_out_seq[ori_seq];
                                    m_C_to_D[atgc.find(sub_seq.strsub[i,1])]+=m_out-Seq[ori_seq];
                                }
                            }  
                        }
                    } else {
                        indel='WT';
                        m_count_list[0]+=m_out_seq[ori_seq];
                    }
                }
            String m_out_seq_read= (String) m_out_seq[ori_seq]; //type casting
            m_emboss_dict[ori_seq]=[sub_wt,sub_sym,sub_seq,m_out_seq_read,indel]; // map String to array
            emboss_wt='';
            emboss_seq='';
            }
        }
    }
    fclose(fp);
    return;
}  

struct less_second {
    typedef pair<string, int> type;
    bool operator () (type const& a, type const &b) const {
        return a.second < b.second;
    }
};
 
vector<pair<string, int> > sorted_map (Sequence m_out_seq) {
    vector<pair<string, int> > cp(m_out_seq.begin(), m_out_seq.end());
    sort(cp.begin(), cp.end(), less_second());
    return cp;
}

void sort (m_out_seq) {
    typedef map <String, int>::iterarot it; // sorting
    String m_sorted_list[m_out_seq.size()];
    int i=0;
    for (it = m_out_seq.begin(); it != myMap.end() ; it++) {
        m_sorted_list[i]=(*it)->first; //??
    }
    return;
}
    
void write_emboss_dict (m_sorted_list, m_emboss_dict) {
    int i;
    FILE* pf;

    pf=fopen("dict.txt",'w');
    if ( pf == NULL) {
        printf("file write error in emboss_dict");
    } else {
        for (i=0, i<m_emboss_dict.size(), i++) {
            fprintf( pf, "%s\t%s", m_sorted_list[i], m_emboss_dict[m_sorted_list[i]]);
        }
    }
    fclose(pf);
    return;
}

void write_emboss_result (m_sorted_list, m_emboss_dict) {
    int i;
    String order;
    FILE* pf;

    pf=fopen("EMBOSS_result.txt",'w');
    if (pf == NULL) {
        printf("file write error in EMBOSS_result.txt");
    } else {
        for (i=0, i<m_emboss_dict.size(), i++) {
            order=m_sorted_list[i];
            fprintf(pf, "%s\t%s\t%s\n%s\n%s\n\n", m_emboss_dict[order][0], m_emboss_dict[order][3], m_emboss_dict[order][4],m_emboss_dict[order][1],m_emboss_dict[order][2]);
        }
    fclose(pf);
    return;
} 
 
void write_count(m_all_count, m_primer_count, m_filted_count, m_count_list) {
    FILE* pf;
    
    pf=fopen("Count.txt",'w');
    if (pf == NULL) {
        printf("file write error in count.txt");
    } else {
        fprintf("all count\t%d\nprimer count\t%d\nfilted count\t%d\n",m_all_count,m_primer_count,m_primer_count-m_filted_count);
        fprintf("WT\t%d\nTotal Substitutions\t%d\nC to D substitutions\t%d\nInsertions\t%d\nDeletions\t%d\n", m_count_list[0], m_count_list[1], m_count_list[2], m_count_list[3], m_count_list[4]);
    }
    fclose(pf);
    return;
}

void write_total_substitution (m_seq_range, target_seq, m_substitution_position, m_substitution_pattern) {
    int i=0, c=0;
    FILE* fp;
    
    pf=fopen("Total_Substitutions.txt",'w');
    if (pf==NULL) {
        printf("file write error in Total_Substitutions.txt");
    } else {
        fprintf("Site\tWT_seq\tSubstitution_count\tA\tT\tG\tC\n");
        c=0;
        for (i=(m_seq_range.find(target_seq)+17)*(-1), i<(m_seq_range.find(target_seq)+17)*(-1)+m_seq_range.size()+1) {
            if (i==0) {
                continue;
            } else {
                c+=1
                fprintf("%d\t%s\t%d\t%d\t%d\t%d\t%d\n", i, m_seq_range.strsub(c-1,1),m_substitution_pattern[c-1],m_substitution_position[0][c-1],m_substitution_position[1][c-1],m_substitution_position[2][c-1],m_substitution_position[3][c-1]);
            }
        }
    }
    fclose(pf);
    return;
}

void write_C_to_D_substitution(m_seq_range, target_seq, m_C_to_D) {
    int i=0, c=0;
    FILE* fp;
    pf=fopne("C_to_D_Substitutions.xlsx", 'w');
    if (pf==NULL) {
        printf("file write error in C_to_D_Substitutions.txt");
    } else {
        fprintf("Site\tWT_seq\tSubstitution_count\tA\tT\tG\tC\n");
        c=0;
        for (i=(m_seq_range.find(target_seq)+17)*(-1), i<(m_seq_range.find(target_seq)+17)*(-1)+m_seq_range.size()+1) {
            if (i==0) {
                continue;
            } else {
                c+=1
                fprintf("%d\t%s\t%d\t%d\t%d\t%d\t%d\n", i, m_seq_range.strsub(c-1,1),m_C_to_D[0][c-1]+m_C_to_D[1][c-1]+m_C_to_D[2][c-1]+m_C_to_D[3][c-1],m_C_to_D[0][c-1],m_C_to_D[1][c-1],m_C_to_D[2][c-1],m_C_to_D[3][c-1]);
            }
        }
    }
    fclose(pf);
    return;
}

