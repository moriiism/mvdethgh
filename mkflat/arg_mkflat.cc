#include "arg_mkflat.h"

// public

void ArgValMkflat::Init(int argc, char* argv[])
{
    progname_ = "mkflat";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 7;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    data_list_      = argv[iarg]; iarg++;
    subimg_dat_     = argv[iarg]; iarg++;
    time_dat_       = argv[iarg]; iarg++;
    vel_dat_        = argv[iarg]; iarg++;
    res_dat_        = argv[iarg]; iarg++;
    outdir_         = argv[iarg]; iarg++;
    outfile_head_   = argv[iarg]; iarg++;
}

void ArgValMkflat::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__, g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n", __func__, progname_.c_str());
    fprintf(fp, "%s: data_list_      : %s\n", __func__, data_list_.c_str());
    fprintf(fp, "%s: subimg_dat_     : %s\n", __func__, subimg_dat_.c_str());
    fprintf(fp, "%s: time_dat_       : %s\n", __func__, time_dat_.c_str());
    fprintf(fp, "%s: vel_dat_        : %s\n", __func__, vel_dat_.c_str());
    fprintf(fp, "%s: res_dat_        : %s\n", __func__, res_dat_.c_str());
    fprintf(fp, "%s: outdir_         : %s\n", __func__, outdir_.c_str());
    fprintf(fp, "%s: outfile_head_   : %s\n", __func__, outfile_head_.c_str());
}

// private

void ArgValMkflat::Null()
{
    progname_     = "";
    data_list_    = "";
    subimg_dat_   = "";
    time_dat_     = "";
    vel_dat_      = "";
    res_dat_      = "";
    outdir_       = "";
    outfile_head_ = "";
}

void ArgValMkflat::SetOption(int argc, char* argv[], option* long_options)
{
    if(0 < g_flag_verbose){
        MPrintInfo("start...");
    }
    // option default
    g_flag_debug   = 0;
    g_flag_help    = 0;
    g_flag_verbose = 0;
    while (1) {
        int option_index = 0;
        int retopt = getopt_long(argc, argv, "dhv",
                                 long_options, &option_index);
        if(-1 == retopt)
            break;
        switch (retopt) {
        case 0:
            // long option
            break;
        case 'd':
            g_flag_debug = atoi(optarg);
            printf("%s: g_flag_debug = %d\n", __func__, g_flag_debug);
            break;
        case 'h':
            g_flag_help = atoi(optarg);
            printf("%s: g_flag_help = %d\n", __func__, g_flag_help);
            if(0 != g_flag_help){
                Usage(stdout);
            }
            break;
        case 'v':
            g_flag_verbose = atoi(optarg);
            printf("%s: g_flag_verbose = %d\n", __func__, g_flag_verbose);
            break;
        case '?':
            printf("%s: retopt (= %c) is invalid flag.\n",
                   __func__, retopt);
            Usage(stdout);
            break;
        default:
            printf("%s: error: getopt returned character code 0%o ??\n", __func__, retopt);
            abort();
        }
    }
    if(0 < g_flag_verbose){
        MPrintInfo("done.");
    }
}


void ArgValMkflat::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "data_list  subimg_dat  time_dat  vel_dat  res_dat  outdir  outfile_head  \n",
            progname_.c_str());
    abort();
}
