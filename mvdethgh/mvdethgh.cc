#include "mir_math.h"
#include "mi_str.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mvdethgh.h"
#include "sub_mvdethgh.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    double time_st = MiTime::GetTimeSec();
    
    ArgValMvdethgh* argval = new ArgValMvdethgh;
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


    string* line_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetDatalist(), &line_arr, &nline);
    string* fitsfile_arr = new string[nline];
    double* time_arr     = new double[nline];
    long ntime = nline;
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(line_arr[iline], &nsplit, &split_arr);
        fitsfile_arr[iline] = split_arr[0];
        time_arr[iline] = atof(split_arr[1].c_str());
        MiStr::DelSplit(split_arr);
    }

    // show datalist
    for(int itime = 0; itime < ntime; itime ++){
        printf("%s  %e\n", fitsfile_arr[itime].c_str(), time_arr[itime]);
    }
    
    printf("--- img_info_in ---\n");
    MifImgInfo* img_info_in = new MifImgInfo;
    img_info_in->Load(argval->GetImgInfoDat());
    img_info_in->PrintInfo();
    printf("=== img_info_in ===\n");

    double** data_arr = new double* [ntime];
    int bitpix = 0;
    for(int itime = 0; itime < ntime; itime ++){
        MifFits::InFitsImageD(fitsfile_arr[itime], img_info_in,
                              &bitpix, &data_arr[itime]);
    }

    MifImgInfo* img_info_rec = new MifImgInfo;
    img_info_rec->InitSetCube(1, 1, 1, img_info_in->GetNaxesArrElm(0),
                              img_info_in->GetNaxesArrElm(1), ntime);
    img_info_rec->PrintInfo();
    double* rec_arr = new double [img_info_rec->GetNpixelTotal()];
    for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
        rec_arr[iarr] = 0.0;
    }
    
    // theta : 0<= theta <= pi/2

    // phi   : 0<= phi <= 2 pi
    // psi   : 0<= psi <= 2 pi
    // rho   : 0<= rho <= 50 

    long ntheta = 100;
    long nphi   = 100;
    long npsi   = 100;
    long nrho   = 100;
    double delta_theta = M_PI / 2. / ntheta;
    double delta_phi   = 2 * M_PI / nphi;
    double delta_psi   = 2 * M_PI / npsi;
    double delta_rho   = (double) img_info_in->GetNaxesArrElm(0) / (double) nrho;    

    MifImgInfo* img_info_cube = new MifImgInfo;
    img_info_cube->InitSetCube(1, 1, 1, nrho, nphi, npsi);
    img_info_cube->PrintInfo();
    
    double* cube_arr = new double [nrho * nphi * npsi];
    for(long iarr = 0; iarr < nrho * nphi * npsi; iarr ++){
        cube_arr[iarr] = 0.0;
    }

    double posx_lo = 0.0;
    double posy_lo = 0.0;
    double posz_lo = 0.0;
    long nposx = img_info_in->GetNaxesArrElm(0);
    long nposy = img_info_in->GetNaxesArrElm(1);
    long nposz = ntime;

    printf("nposx = %ld\n", nposx);
    printf("nposy = %ld\n", nposy);
    printf("ntime = %ld\n", ntime);

    for(long itheta = 0; itheta < ntheta; itheta ++){
        double theta = delta_theta * itheta;
        for(long iarr = 0; iarr < nrho * nphi * npsi; iarr ++){
            cube_arr[iarr] = 0.0;
        }
        for(long iarr = 0; iarr < img_info_rec->GetNpixelTotal(); iarr++){
            rec_arr[iarr] = 0.0;
        }
        for(long itime = 0; itime < ntime; itime ++){
            for(long iposx = 0; iposx < nposx; iposx ++){
                for(long iposy = 0; iposy < nposy; iposy ++){
                    double xval = iposx;
                    double yval = iposy;
                    double zval = itime;

                    if(data_arr[itime][iposx + iposy * nposx] < 1e-10){
                        continue;
                    }

                    for(long ipsi = 0; ipsi < npsi; ipsi ++){
                        double psi = delta_psi * ipsi;
                        double rho_cos_phi = xval * cos(psi) + yval * sin(psi);
                        double rho_sin_phi =
                            -1 * xval * sin(psi) * cos(theta)
                            + yval * cos(psi) * cos(theta)
                            + zval * sin(theta);
                        double rho = sqrt(rho_cos_phi * rho_cos_phi + rho_sin_phi * rho_sin_phi);
                        double phi = 0.0;
                        if(rho_sin_phi >= 0.0){
                            phi = acos( rho_cos_phi /
                                        sqrt( rho_cos_phi * rho_cos_phi + rho_sin_phi * rho_sin_phi ) );
                        } else {
                            phi = acos( rho_cos_phi /
                                        sqrt( rho_cos_phi * rho_cos_phi + rho_sin_phi * rho_sin_phi ) )
                                + M_PI;
                        }
                        long irho = (long) floor( rho / delta_rho );
                        long iphi = (long) floor( phi / delta_phi );
                        long i_rho_phi_psi = irho + iphi * nrho + ipsi * (nrho * nphi);
                        cube_arr[i_rho_phi_psi] += data_arr[itime][iposx + iposy * nposx];
                    }
                }
            }
        }
        char tag[kLineSize];
        sprintf(tag, "theta_%2.2ld", itheta);
        int naxis = 2;
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_cube->GetNaxesArr(),
                              cube_arr);
        // get max
        double max_arr = 0.0;
        long index = 0;
        for(long iarr = 0; iarr < nrho * nphi * npsi; iarr ++){
            if(cube_arr[iarr] > max_arr){
                max_arr = cube_arr[iarr];
                index = iarr;
            }
        }

        long ipsi = index / (nrho * nphi);
        long index2 = index % (nrho * nphi);
        long iphi = index2 / nrho;
        long irho = index2 % nrho;

        double rho = irho * delta_rho;
        double phi = iphi * delta_phi;
        double psi = ipsi * delta_psi;

        printf("itheta = %ld, theta = %e, rho = %e, phi = %e, psi = %e \n",
               itheta, theta, rho, phi, psi);

        long ntpar = 1000;
        double delta_tpar = 1.0;
        double tpar_lo = -10;
        for(long itpar = 0; itpar < ntpar; itpar ++){
            double tpar = tpar_lo + itpar * delta_tpar;
            double xval = rho * (cos(psi) * cos(phi) - sin(psi) * cos(theta) * sin(phi) )
                + tpar * sin(psi) * sin(theta);
            double yval = rho * (sin(psi) * cos(phi) + cos(psi) * cos(theta) * sin(phi) )
                - tpar * cos(psi) * sin(theta);
            double zval = rho * sin(theta) * sin(phi) + tpar * cos(theta);

            long iposx = (long) xval;
            long iposy = (long) yval;
            long iposz = (long) zval;

            // printf("itpar = %ld\n", itpar);
            long index = iposx + iposy * nposx + iposz * (nposx * nposy);
            // printf("index = %ld, total = %ld\n", index, nposx * nposy * nposz);

            if(iposx < 0 || nposx - 1 < iposx){
                continue;
            }
            if(iposy < 0 || nposy - 1 < iposy){
                continue;
            }
            if(iposz < 0 || nposz - 1 < iposz){
                continue;
            }
            printf("fill: index = %ld\n", index);
            rec_arr[index] = 1.0;
        }

        sprintf(tag, "rec_%2.2ld", itheta);
        MifFits::OutFitsCubeD(argval->GetOutdir(), argval->GetOutfileHead(), tag,
                              3, bitpix,
                              img_info_rec->GetNaxesArr(),
                              rec_arr);
        
        // get (rho, phi, psi)

        // long ipsi = index / (nrho * nphi)
        
        
    }


    
    return status_prog;
}
