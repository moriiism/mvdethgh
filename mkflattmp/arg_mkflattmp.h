#ifndef MORIIISM_MVDETHGH_MKFLATTMP_ARG_MKFLATTMP_H_
#define MORIIISM_MVDETHGH_MKFLATTMP_ARG_MKFLATTMP_H_

#include "mi_base.h"

class ArgValMkflattmp : public MiArgBase{
public:
    ArgValMkflattmp() :
        MiArgBase(),
        progname_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkflattmp(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_MKFLATTMP_ARG_MKFLATTMP_H_
