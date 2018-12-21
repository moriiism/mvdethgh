#include "mir_math.h"
#include "mi_str.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_simmv.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValSimmv* argval = new ArgValSimmv;
    argval->Init(argc, argv);
    argval->Print(stdout);

    char logfile[kLineSize];
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        system(cmd);
    }
    sprintf(logfile, "%s/%s_%s.log",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str(),
            argval->GetProgname().c_str());
    FILE* fp_log = fopen(logfile, "w");
    MiIolib::Printf2(fp_log, "-----------------------------\n");
    argval->Print(fp_log);

    
    printf("--- img_info_in ---\n");
    MifImgInfo* img_info_in = new MifImgInfo;
    img_info_in->Load(argval->GetImgInfoDat());
    img_info_in->PrintInfo();
    printf("=== img_info_in ===\n");

    int nimg = argval->GetNimg();
    double** sim_arr = new double* [nimg];
    for(int iimg = 0; iimg < nimg; iimg ++){
        sim_arr[iimg] = new double [img_info_in->GetNpixelTotal()];
        for(long iarr = 0; iarr < img_info_in->GetNpixelTotal(); iarr ++){
            sim_arr[iimg][iarr] = 0.0;
        }
    }

    for(int iimg = 0; iimg < nimg; iimg ++){

        double xvel = 1.0;
        double yvel = -0.5;
        double xpos_st = 20.0;
        double ypos_st = 20.0;
        double xpos = xpos_st + xvel * iimg;
        double ypos = ypos_st + yvel * iimg;

        long iposx = (long) floor(xpos);
        long iposy = (long) floor(ypos);
        long iarr = iposx + iposy * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;

        iarr = (iposx + 1) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
        iarr = (iposx + 1) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
        iarr = (iposx + 1) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;

        iarr = (iposx + 0) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
        iarr = (iposx + 0) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
        iarr = (iposx + 0) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
        
        iarr = (iposx - 1) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
        iarr = (iposx - 1) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
        iarr = (iposx - 1) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
        sim_arr[iimg][iarr] += 10;
    }
//
//    for(int iimg = 0; iimg < nimg; iimg ++){
//
//        double xvel = -1.0;
//        double yvel = 2.0;
//        double xpos_st = 20.0;
//        double ypos_st = 20.0;
//        double xpos = xpos_st + xvel * iimg;
//        double ypos = ypos_st + yvel * iimg;
//
//        long iposx = (long) floor(xpos);
//        long iposy = (long) floor(ypos);
//        long iarr = iposx + iposy * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//
//        iarr = (iposx + 1) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//        iarr = (iposx + 1) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//        iarr = (iposx + 1) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//
//        iarr = (iposx + 0) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//        iarr = (iposx + 0) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//        iarr = (iposx + 0) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//        
//        iarr = (iposx - 1) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//        iarr = (iposx - 1) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//        iarr = (iposx - 1) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 5.0;
//    }
    
//    for(int iimg = 0; iimg < nimg; iimg ++){
//
//        double xvel = 0.0;
//        double yvel = 0.0;
//        double xpos_st = 20.0;
//        double ypos_st = 20.0;
//        double xpos = xpos_st + xvel * iimg;
//        double ypos = ypos_st + yvel * iimg;
//
//        long iposx = (long) floor(xpos);
//        long iposy = (long) floor(ypos);
//        long iarr = iposx + iposy * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//
//        iarr = (iposx + 1) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//        iarr = (iposx + 1) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//        iarr = (iposx + 1) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//
//        iarr = (iposx + 0) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//        iarr = (iposx + 0) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//        iarr = (iposx + 0) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//        
//        iarr = (iposx - 1) + (iposy + 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//        iarr = (iposx - 1) + (iposy + 0) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//        iarr = (iposx - 1) + (iposy - 1) * img_info_in->GetNaxesArrElm(0);
//        sim_arr[iimg][iarr] += 1.0;
//    }
    
    
    for(int iimg = 0; iimg < nimg; iimg ++){
        char tag[kLineSize];
        sprintf(tag, "sim_%2.2d", iimg);
        int bitpix = -32;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                               2, bitpix,
                               img_info_in->GetNaxesArr(),
                               sim_arr[iimg]);
    }
    
    return status_prog;
}
