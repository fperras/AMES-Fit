#include "powder_calc.hpp"
using namespace std;

void reverse_axis(int TD, vector<float> &bin, vector<float> &ppm){
    //Function used to change the order of the chemical shift axis when reading data to conform with the order used in AMES-Fit.
    int i,j=TD-1;
    float* binI;
    float* ppmI;
    binI = (float*)malloc(TD * sizeof(float));
    ppmI = (float*)malloc(TD * sizeof(float));

    for(i=0;i<TD;i++){
        ppmI[i]=ppm[j];
        binI[i]=bin[j];
        j--;
    }
    for(i=0;i<TD;i++){
        ppm[i]=ppmI[i];
        bin[i]=binI[i];
    }
    free(binI);
    free(ppmI);
}

float load_dat(char* filename, int TD, vector<float> &bin, vector<float> &ppm){
    //Function used to load an experimental spectrum to memory.
    //Intensities are stored under bin while the ppm values of these points under ppm.
    //A specific TD can be specified to downscale a spectrum and improve efficiency.
    //This is done by averaging all points that fall within a given window.
    //Data must be in the SOLIDS format
    FILE *fp;
    fp=fopen(filename,"r");
    char buffer[256];
    float high, low, vs;
    int TDI, i, j;

    for (i=0; i<TD; i++){
        bin[i]=0.;
    }

    //Reading the header to the SOLIDS file
    fgets(buffer, 128, fp);
    sscanf(buffer," %d",&TDI);

    fgets(buffer, 128, fp);
    sscanf(buffer," %f",&vs);

    fgets(buffer, 128, fp);
    fgets(buffer, 128, fp);
    sscanf(buffer," %f",&high);

    fgets(buffer, 128, fp);
    sscanf(buffer," %f",&low);

    float* binI;
    float* ppmI;
    binI = (float*)malloc(TDI * sizeof(float));
    ppmI = (float*)malloc(TDI * sizeof(float));

    //Creating the shift axis
    float ppm_incrI=(high-low)/((TDI+1)*vs);
    ppmI[0]=high/vs-ppm_incrI/2.;
    for(i=1;i<TDI;i++){
        ppmI[i]=ppmI[i-1]-ppm_incrI;
    }

    //Reading the intensities
    for(i=0;i<TDI;i++){
        fgets(buffer, 128, fp);
        sscanf(buffer," %f",&binI[i]);
    }

    float ppm_incr=(high-low)/((TD+1)*vs);
    ppm[0]=high/vs-0.5*ppm_incr;
    vs=vs*1.e6;
    for(i=1;i<TD;i++){
        ppm[i]=ppm[i-1]-ppm_incr;
    }

    //Scaling back the intensities to avoid overflow with large downscaling ratios
    bin[0]=binI[0]*(TDI/TD);
    bin[TD-1]=binI[TDI-1]*(TDI/TD);
    int counter[TD];

    //downscaling
    for(i=1;i<TD-1;i++){
        counter[i]=bin[i]=0;
        for(j=0;j<TDI-1;j++){
            if((ppm[i]<=ppmI[j])&&(ppm[i-1]>=ppmI[j])){
                counter[i]++;
                bin[i]=bin[i]+binI[j];
            }
        }
        bin[i]=bin[i]/counter[i];
    }

    //Setting the maximum intensity to 1.0
    normalize_spectrum(bin,TD);
    reverse_axis(TD,bin,ppm);

    free(ppmI);
    free(binI);
    fclose(fp);
    return vs;
}


