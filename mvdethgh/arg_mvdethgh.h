#ifndef MORIIISM_MVDETHGH_ARG_MVDETHGH_H_
#define MORIIISM_MVDETHGH_ARG_MVDETHGH_H_

#include "mi_base.h"

class ArgValMvdethgh : public MiArgBase{
public:
    ArgValMvdethgh() :
        MiArgBase(),
        progname_(""),
        datalist_(""),
        img_info_dat_(""),
        par_dat_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMvdethgh(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDatalist() const {return datalist_;};
    string GetImgInfoDat() const {return img_info_dat_;};
    string GetParDat() const {return par_dat_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string datalist_;
    string img_info_dat_;
    string par_dat_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_ARG_MVDETHGH_H_
