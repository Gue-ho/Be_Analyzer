#include <zlib.h>
#include <sys/stat.h>

#include "be_analyzer.h"
#include "be_fastq-lib.h"
#include "emscripten_progress.h"

#define VERSION "1.3"
#define SVNREV 1

int BeAnalyzer::be_fastq_join(const char* fastq_file1, const char* fastq_file2) {
    int debug=0;
    char c;
    int mismatch = 0;
    char *in[3] = {0,0,0};
    char *out[5];
    char *orep=NULL;
    int in_n = 2;
    int threads = 1;				// not really necessary
    char verify='\0';

    int i;
    int mino = 6;
    int pctdiff = 8;				// this number tested well on exome data... tweak for best results
    bool omode = false;
    char *bfil = NULL;
    bool norevcomp = false;
    bool allow_ex = false;

    in[0] = (char *)fastq_file1;
    in[1] = (char *)fastq_file2;

    report_progress(0, 100, 0, 100, "Loading files...");

    gzFile fin[2];
    for (i = 0; i < in_n; ++i) {
        FILE* fp = fopen(in[i], "r");
        setvbuf(fp, NULL, _IOFBF, 4194304);
        fin[i] = gzdopen(fileno(fp), "r");
        if (!fin[i]) {
            fprintf(stderr, "Error opening file '%s': %s\n",in[i], strerror(errno));
            return 1;
        }
    }

    // some basic validation of the file formats
    {
        for (i=0;i<in_n;++i) {
            char c=gzgetc(fin[i]);
            if (c != '@')  {
                fprintf(stderr, "%s doesn't appear to be a fastq file (%c)\n", in[i], c);
                return 1;
            }
            gzungetc(c, fin[i]);
        }
    }

    struct fq fq[2];
    meminit(fq);

    int nrec=0;
    int nerr=0;
    int nok=0;
    int joincnt=0;
    double tlen=0;
    double tlensq=0;
    int read_ok;

    struct fq rc;
    meminit(rc);

    struct stat st;
    stat(in[0], &st);
    double total_size = (double)st.st_size;
    //double prev_pos = -1;

    // read in 1 record from each file
    while ((read_ok=read_fq(fin[0], nrec, &fq[0]))) {
        double current_pos = (double)gzoffset(fin[0]);
        if (nrec % 1000 == 0) {
            report_progress(current_pos, total_size, 0, 30, "Running fastq-join...");
            //prev_pos = current_pos;
        }
        for (i=1;i<in_n;++i) {
            int mate_ok=read_fq(fin[i], nrec, &fq[i]);
            if (read_ok != mate_ok) {
                fprintf(stderr, "# of rows in mate file '%s' doesn't match primary file, quitting!\n", in[i]);
                return 1;
            }
        }

        ++nrec;
        if (read_ok < 0) continue;

        if (!norevcomp) {
            revcomp(&rc, &fq[1]);
        }

        int maxo = min(fq[0].seq.n, rc.seq.n);
        int bestscore=INT_MAX;
        int besto=-1;
        for (i=mino; i <= maxo; ++i) {
            int mind = (pctdiff * i) / 100;
            int d;
            d=hd(fq[0].seq.s+fq[0].seq.n-i, rc.seq.s, i);
            if (d <= mind) {
                // squared-distance over length, probably can be proven better (like pearson's)
                int score = (1000*(d*d+1))/i;
                if (score < bestscore) {
                    bestscore=score;
                    besto=i;
                }
            }
        }

        int hasex=0;

        int olen = besto-hasex;

        if (besto > 0) {
            ++joincnt;

            tlen+=olen;
            tlensq+=olen*olen;

            char *sav_fqs=NULL, *sav_rcs;
            char *sav_fqq, *sav_rcq;

            for (i = 0; i < besto; ++i ) {
                int li = fq[0].seq.n-besto+i;
                int ri = i;
                if (fq[0].seq.s[li] == rc.seq.s[ri]) {
                    fq[0].qual.s[li] = max(fq[0].qual.s[li], rc.qual.s[ri]);
                } else {
                    if (fq[0].qual.s[li] > rc.qual.s[ri]) {
                        fq[0].qual.s[li] = 33+min(fq[0].qual.s[li],max(fq[0].qual.s[li]-rc.qual.s[ri],3));
                    } else {
                        fq[0].seq.s[li] = rc.seq.s[ri];
                        fq[0].qual.s[li] = 33+min(rc.qual.s[ri],max(rc.qual.s[ri]-fq[0].qual.s[li],3));
                    }
                }
            }
            update_seq(string(fq[0].seq.s) + string(rc.seq.s+besto));
        }
    }

    /*
    double dev = sqrt((((double)joincnt)*tlensq-pow((double)tlen,2)) / ((double)joincnt*((double)joincnt-1)) );
    printf("Total reads: %d\n", nrec);
    printf("Total joined: %d\n", joincnt);
    printf("Average join len: %.2f\n", (double) tlen / (double) joincnt);
    printf("Stdev join len: %.2f\n", dev);
    printf("Version: %s.%d\n", VERSION, SVNREV);
    */
    return 0;
}
