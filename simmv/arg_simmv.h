#ifndef MORIIISM_MVDETHGH_ARG_SIMMV_H_
#define MORIIISM_MVDETHGH_ARG_SIMMV_H_

#include "mi_base.h"

class ArgValSimmv : public MiArgBase{
public:
    ArgValSimmv() :
        MiArgBase(),
        progname_(""),
        img_info_dat_(""),
        nimg_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValSimmv(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetImgInfoDat() const {return img_info_dat_;};
    int    GetNimg() const {return nimg_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string img_info_dat_;
    int    nimg_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_ARG_SIMMV_H_
