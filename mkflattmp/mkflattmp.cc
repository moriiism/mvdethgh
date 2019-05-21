#include "mi_time.h"
#include "mvdethghlib.h"
#include "arg_mkflattmp.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMkflattmp* argval = new ArgValMkflattmp;
    argval->Init(argc, argv);
    argval->Print(stdout);
    
    FILE* fp_log = NULL;
    OpenLogfile(argval->GetOutdir(),
                argval->GetOutfileHead(),
                argval->GetProgname(),
                &fp_log);
    argval->Print(fp_log);

    MirRootTool* root_tool = new MirRootTool;
    root_tool->InitTCanvas("pub");

    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, 100, 100);
    
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(100, 0.0, 1.0, 100, 0.0, 1.0);
    
    double val_lo = 0.0;
    double delta_bin = 0.0001;
    for(long ibin = 0; ibin < 10000; ibin ++){
        double val = val_lo + (ibin + 0.5) * delta_bin;
        for(long ibinx = 0; ibinx < hd2d->GetHi2d()->GetNbinX(); ibinx ++){
            double xval = hd2d->GetHi2d()->GetBinCenterXFromIbinX(ibinx);
            double yval = val / xval;
            printf("xval, yval = %e, %e\n", xval, yval);
            if(hd2d->GetHi2d()->GetLoX() < xval && xval < hd2d->GetHi2d()->GetUpX() &&
               hd2d->GetHi2d()->GetLoY() < yval && yval < hd2d->GetHi2d()->GetUpY() ){
                hd2d->Fill(xval, yval);
            }
        }
    }

    hd2d->MkTH2Fig("temp.png", root_tool);
    
    int bitpix = -32;
    MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(), "tmp",
                           2, bitpix,
                           img_info->GetNaxesArr(),
                           hd2d->GetOvalArr()->GetVal());

    double time_ed = MiTime::GetTimeSec();
    printf("time_ed - time_st = %e\n", time_ed - time_st);
    
    return status_prog;
}
