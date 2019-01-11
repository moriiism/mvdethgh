#ifndef MORIIISM_MVDETHGH_MKPSF_ARG_MKPSF_H_
#define MORIIISM_MVDETHGH_MKPSF_ARG_MKPSF_H_

#include "mi_base.h"

class ArgValMkpsf : public MiArgBase{
public:
    ArgValMkpsf() :
        MiArgBase(),
        progname_(""),
        data_list_(""),
        subimg_dat_(""),
        func_(""),
        par_file_(""),
        nbin_kernel_half_(0),
        nbin_psf_half_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkpsf(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataList() const {return data_list_;};
    string GetSubimgDat() const {return subimg_dat_;};
    string GetFunc() const {return func_;};
    string GetParFile() const {return par_file_;};
    int    GetNbinKernelHalf() const {return nbin_kernel_half_;};
    int    GetNbinPsfHalf() const {return nbin_psf_half_;};    
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string data_list_;
    string subimg_dat_;
    string func_;
    string par_file_;
    int    nbin_kernel_half_;
    int    nbin_psf_half_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_MKPSF_ARG_MKPSF_H_
