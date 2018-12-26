#ifndef MORIIISM_MVDETHGH_ARG_MVDETHGH_H_
#define MORIIISM_MVDETHGH_ARG_MVDETHGH_H_

#include "mi_base.h"

class ArgValMvdethgh : public MiArgBase{
public:
    ArgValMvdethgh() :
        MiArgBase(),
        progname_(""),
        data_list_(""),
        subimg_dat_(""),
        time_dat_(""),
        vel_dat_(""),
        res_dat_(""),
        sig_(0.0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMvdethgh(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataList() const {return data_list_;};
    string GetSubimgDat() const {return subimg_dat_;};
    string GetTimeDat() const {return time_dat_;};
    string GetVelDat() const {return vel_dat_;};
    string GetResDat() const {return res_dat_;};
    double GetSig() const {return sig_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string data_list_;
    string subimg_dat_;
    string time_dat_;
    string vel_dat_;
    string res_dat_;
    double sig_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_ARG_MVDETHGH_H_
