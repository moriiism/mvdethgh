#ifndef MORIIISM_MVDETHGH_MVDETSUM_ARG_MVDETSUM_H_
#define MORIIISM_MVDETHGH_MVDETSUM_ARG_MVDETSUM_H_

#include "mi_base.h"

class ArgValMvdetsum : public MiArgBase{
public:
    ArgValMvdetsum() :
        MiArgBase(),
        progname_(""),
        data_list_(""),
        subimg_dat_(""),
        time_dat_(""),
        vel_dat_(""),
        npsi_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMvdetsum(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataList() const {return data_list_;};
    string GetSubimgDat() const {return subimg_dat_;};
    string GetTimeDat() const {return time_dat_;};
    string GetVelDat() const {return vel_dat_;};
    int    GetNpsi() const {return npsi_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string data_list_;
    string subimg_dat_;
    string time_dat_;
    string vel_dat_;
    int    npsi_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_MVDETHGH_MVDETSUM_ARG_MVDETSUM_H_
