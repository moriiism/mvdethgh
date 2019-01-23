#ifndef MORIIISM_MVDETHGH_MKDETIMG_ARG_MKDETIMG_H_
#define MORIIISM_MVDETHGH_MKDETIMG_ARG_MKDETIMG_H_

#include "mi_base.h"

class ArgValMkdetimg : public MiArgBase{
public:
    ArgValMkdetimg() :
        MiArgBase(),
        progname_(""),
        data_list_(""),
        par_dat_(""),
        nbin_detimg_half_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkdetimg(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataList() const {return data_list_;};
    string GetParDat() const {return par_dat_;};
    int    GetNbinDetimgHalf() const {return nbin_detimg_half_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string data_list_;
    string par_dat_;
    int    nbin_detimg_half_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_MKDETIMG_ARG_MKDETIMG_H_
