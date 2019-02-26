#ifndef MORIIISM_MVDETHGH_MVDETHGHXYT_ARG_MVDETHGHXYT_H_
#define MORIIISM_MVDETHGH_MVDETHGHXYT_ARG_MVDETHGHXYT_H_

#include "mi_base.h"

class ArgValMvdethghxyt : public MiArgBase{
public:
    ArgValMvdethghxyt() :
        MiArgBase(),
        progname_(""),
        data_list_(""),
        vel_dat_(""),
        res_dat_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMvdethghxyt(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataList() const {return data_list_;};
    string GetVelDat() const {return vel_dat_;};
    string GetResDat() const {return res_dat_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string data_list_;
    string vel_dat_;
    string res_dat_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_MVDETHGHXYT_ARG_MVDETHGHXYT_H_
