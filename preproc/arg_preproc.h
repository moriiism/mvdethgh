#ifndef MORIIISM_MVDETHGH_PREPROC_ARG_PREPROC_H_
#define MORIIISM_MVDETHGH_PREPROC_ARG_PREPROC_H_

#include "mi_base.h"

class ArgValPreproc : public MiArgBase{
public:
    ArgValPreproc() :
        MiArgBase(),
        progname_(""),
        data_list_(""),
        subimg_dat_(""),
        psf_dat_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValPreproc(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataList() const {return data_list_;};
    string GetSubimgDat() const {return subimg_dat_;};
    string GetPsfDat() const {return psf_dat_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string data_list_;
    string subimg_dat_;
    string psf_dat_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_PREPROC_ARG_PREPROC_H_
