#include <string>
#include <algorithm>

using namespace std;

class BeAnalyzer {
    private:
        BeAnalyzer();
        void seq_standard(String, String, int, int);
        String m_s_seq, m_seq_range;
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

    for (i=22; i<wt_seq.size()-23; i++) {
        if (wt_seq.substr(i, 23).compare(target_seq) == 0) {
            start_pos = i+ 17 - end_range;
            end_pos = i + 17 + end_range;
            
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
    m_seq_range.substr();
}
