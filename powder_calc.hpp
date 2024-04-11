#include <time.h>
#include <math.h>
#include <vector>
#include "FastExp.h"
#include <gsl/gsl_multimin.h>
using namespace std;

struct par{
    int nspec;
    int fast;
    vector<int> MAS;
    vector<int> sat;
    vector<int> TD;
    int sites;
    vector< vector<float> > bin;
    vector< vector<float> > ppm;
    float S;
    vector<float> vs;
    vector<float> width;
};

int floatcomp(const void* elem1, const void* elem2) {
    if (*(const float*)elem1 < *(const float*)elem2) return -1;
    return *(const float*)elem1 > *(const float*)elem2;
}

float vQ_MAS(float diso, float chi, float eta, float vs, float S, float ct, float c2p) {
    //This function returns the MAS resonance frequency for a specific isochromat
    float hc2p = eta * c2p;
    float F0 = 21. / 16. - 7. / 8. * hc2p + 7. / 48. * hc2p * hc2p;
    float F2 = -9. / 8. + eta * eta / 12. + hc2p - 7. / 24. * hc2p * hc2p;
    float F4 = 5. / 16. - 1. / 8. * hc2p + 7. / 48. * hc2p * hc2p;
    return diso - 9.e6 * chi * chi / (6. * vs * vs * 4. * S * S * (2. * S - 1) * (2. * S - 1)) * (S * (S + 1.) - 0.75) * (F0 * pow(ct, 4.) + F2 * ct * ct + F4);
}

float vQ_stat(float diso, float chi, float eta, float vs, float S, float ct, float c2p) {
     //This function returns the static resonance frequency for a specific isochromat, including second-order quadrupolar coupling and isotropic shift
    float hc2p = eta * c2p;
    float A = -27. / 8. + 9. / 4. * hc2p - 3. / 8. * hc2p * hc2p;
    float B = 30. / 8. - 0.5 * eta * eta - 2. * hc2p + 3. / 4. * hc2p * hc2p;
    float C = -3. / 8. + eta * eta / 3. - 0.25 * hc2p - 3. / 8. * hc2p * hc2p;
    return diso - 9.e6 * chi * chi / (6. * vs * vs * (2. * S * (2. * S - 1.)) * (2. * S * (2. * S - 1.))) * (S * (S + 1.) - 0.75) * (A * pow(ct, 4.) + B * ct * ct + C);
}

float vQ_sat(float diso, float chi, float eta, float vs, float S, float m, float ct, float st, float c2p){
    //returns the first order quadrupolar shift from the transition from m and m-1
    float vQ = (1.-2.*m)*(3.*chi)/(4.*S*(2.*S-1.)) * (0.5*(3.*ct*ct-1) + 0.5*eta*st*st*c2p);
    return diso + vQ*1.e6/vs;
}

float vCSA_stat(float d11, float d22, float d33, float ct, float st, float cp, float sp, float ca, float sa, float cb, float sb, float cg, float sg) {
    //This function returns the resonance frequency contribution to an isochromat from the CSA
    //variables cs and sa are the sine and cosine of the Euler angle alpha, etc.
    float sbsp = sb * sg * ct + ((sa * cb * sg - ca * cg) * st * sp + (-ca * cb * sg - cg * sa) * st * cp);
    float sbcp = sb * cg * ct + ((ca * sg + cg * sa * cb) * st * sp + (sa * sg - ca * cg * cb) * st * cp);
    float cb2 = cb * ct + -sb * sa * st * sp + sb * ca * st * cp;
    return d11 * sbcp * sbcp + d22 * sbsp * sbsp + d33 * cb2 * cb2 - (d11 + d22 + d33) / 3.;
}

void normalize_spectrum(vector<float>& bin, int TD) {
    //Sets the maximum intensity of the spectrum to 1.0
    float max = bin[0];
    int i;
    for (i = 1; i < TD; i++) {
        if (bin[i] > max) max = bin[i];
    }
    for (i = 0; i < TD; i++) {
        bin[i] = bin[i] / max;
    }
}

void Broadening(float LB, float GB, vector<float>& bin, float width, int TD) {
    //Convolutes a spectrum with Gaussian and Laurentzian Broadening
    float* temp_bin = (float*)malloc(TD * sizeof(float));
    float* scaling = (float*)malloc(TD * sizeof(float));
    int i, j;
    float offset;

    for (i = 0; i < TD; i++) {
        temp_bin[i] = 0.;
    }

    if (LB > 0) {
        for (i = 0; i < TD; i++) {
            offset = width * i / TD;
            scaling[i] = 1. / (1. + (4 * offset * offset / (LB * LB)));
        }
        for (i = 0; i < TD; i++) {
            for (j = 0; j < TD; j++) {
                offset = width * abs(i - j) / TD;
                temp_bin[j] = temp_bin[j] + bin[i] * scaling[abs(i - j)];
            }
        }

        for (i = 0; i < TD; i++) {
            bin[i] = temp_bin[i];
        }
    }

    if (GB > 0) {
        for (i = 0; i < TD; i++) {
            offset = width * i / TD;
            scaling[i] = fastExp(-2.77258 * offset * offset / (GB * GB));
        }
        for (i = 0; i < TD; i++) {
            for (j = 0; j < TD; j++) {
                temp_bin[j] = temp_bin[j] + bin[i] * scaling[abs(i - j)];
            }
        }

        for (i = 0; i < TD; i++) {
            bin[i] = temp_bin[i];
        }
    }
    free(temp_bin);

    normalize_spectrum(bin, TD);
}

void simulate_MAS(int TD, vector<float>& bin, vector<float>& ppm, float diso, float chi, float eta, float S, float vs, int fast) {
    //Function to calculate an MAS spectrum

    //parameter declarations
    int N = 32, i, j, k, v, quadrant;
    float l, m, n, ct, st, cp, sp, R, f[3];
    float spacing = ppm[1] - ppm[0];
    int Np = N + 1;
    float* freq;
    freq = (float*)malloc(4 * Np * Np * sizeof(float));

    //initializing the spectrum intensities
    for (v = 0; v < TD; v++) {
        bin[v] = 0.;
    }

    //Looping over the four quadrants of the semi-octahedron for the ASG interpolation.
    //Frequencies are stored under the parameter freq according to the indices of the
    //specific triangle intersections, and the quadrant (face) they are on.
    for (j = 0; j <= N; j++) {
        for (i = 0; (i + j) <= N; i++) {
            for (quadrant = 0; quadrant < 4; quadrant++) {
                k = N - abs(i) - abs(j);
                l = abs(i) / sqrt(i * i + j * j + k * k);
                m = abs(j) / sqrt(i * i + j * j + k * k);
                n = k / sqrt(i * i + j * j + k * k);

                //pre-computing the trigonometric operations
                ct = n;
                st = sqrt(1. - n * n);

                switch (quadrant) {
                    case 0:
                        cp = -l / st;
                        sp = m / st;
                        break;

                    case 1:
                        cp = -l / st;
                        sp = -m / st;
                        break;

                    case 2:
                        cp = l / st;
                        sp = m / st;
                        break;

                    case 3:
                        cp = l / st;
                        sp = -m / st;
                        break;
                }

                if (cp < -1.) {
                    cp = -1;
                }
                if (cp > 1.) {
                    cp = 1;
                }
                if (ct == 1.) {
                    cp = 0.;
                    st = 0.;
                }
                if (ct == 0.) {
                    st = 1.;
                }
                if (cp == 1.) {
                    sp = 0.;
                }
                float c2p = cp * cp - sp * sp;

                //calculating the frequency for the orientation
                freq[quadrant * Np * Np + abs(i) * Np + abs(j)] = vQ_MAS(diso, chi, eta, vs, S, ct, c2p);
            }
        }
    }
    //adding the triangle interpolated "tents" to the spectrum
    int a;
    for (j = 1; j <= N; j++) {
        for (i = 0; i <= (N - abs(j)); i++) {
            int jump;
            for (quadrant = 0; quadrant < 4; quadrant++) {
                for (a = 0; a < 2; a++) {
                    jump = 1;
                    if (a == 0) {
                        f[0] = freq[quadrant * Np * Np + Np * i + j];
                        f[1] = freq[quadrant * Np * Np + Np * i + j - 1];
                        f[2] = freq[quadrant * Np * Np + Np * (i + 1) + j - 1];
                    }

                    else if (a == 1 && (i + j) < (N)) {
                        f[0] = freq[quadrant * Np * Np + Np * i + abs(j)];
                        f[1] = freq[quadrant * Np * Np + Np * (i + 1) + j - 1];
                        f[2] = freq[quadrant * Np * Np + Np * (i + 1) + j];
                    } else
                        jump = 0;

                    k = N - (i) - (j);
                    R = sqrt(i * i + j * j + k * k) / N;
                    qsort(f, 3, sizeof(float), floatcomp);  // in ascending order

                    if (jump) {
                        int v1 = (int)((f[0] - ppm[0]) / spacing);
                        int v3 = (int)((f[2] - ppm[0]) / spacing);
                        int v2 = (int)((f[1] - ppm[0]) / spacing);
                        int vstart = v1;
                        int vend = v3;
                        int vmid = v2;
                        if (vstart < 0) {
                            vstart = 0;
                            v1 = -100e6;
                        }
                        if (vend > TD - 1) {
                            vend = TD - 1;
                            v3 = 100e6;
                        }
                        if (vmid < 0) vmid = 0;
                        if (vmid > TD - 1) vmid = TD - 1;

                        //For initial optimizations this approximate interpolation is used
                        if(fast){
                            for (v = vstart; v < vmid; v++) {
                                bin[v] = bin[v] + (v1!=v2)*(v3!=v1)* (float)(v - v1) / (v2 - v1) / (R * R * R * (v3 - v1));
                            }
                            for (v = vmid; v < vend; v++) {
                                bin[v] = bin[v] + (v3!=v2)*(v3!=v1)*(float)(v3 - v) / (v3 - v2) / (R * R * R * (v3 - v1));
                            }
                        }

                        //Later optimizations consider all ASG possibilities, see the original paper
                       else{
                            for (v=vstart; v<TD-1; v++){
                                if((f[0]>ppm[TD-1]))
                                    break;

                                else if((f[2]<ppm[0]))
                                    break;

                                if((ppm[v]>f[0]) && (ppm[v+1]<=f[1])){
                                    bin[v]=bin[v]+1./(R*R*R)*((ppm[v+1]-ppm[v])*(ppm[v]+ppm[v+1]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[0]) && (ppm[v+1]<=f[1])){
                                    bin[v]=bin[v]+1./(R*R*R)*((ppm[v+1]-f[0])*(ppm[v+1]-f[0]))/((f[2]-f[0])*(f[1]-f[0]));
                                }
                                else if((ppm[v]>f[1]) && (ppm[v]<=f[2]) && (ppm[v+1]>f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*((f[2]-ppm[v])*(f[2]-ppm[v]))/((f[2]-f[0])*(f[2]-f[1]));
                                }
                                else if((ppm[v]>f[0]) && (ppm[v]<=f[1]) && (ppm[v+1]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*(((f[1]-ppm[v])*(f[1]+ppm[v]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]))+((ppm[v+1]-f[1])*(2.0*f[2]-ppm[v+1]-f[1]))/((f[2]-f[0])*(f[2]-f[1])));
                                }
                                else if((ppm[v]>f[0]) && (ppm[v]<=f[1]) && (ppm[v+1]>f[2])){
                                bin[v]=bin[v]+1./(R*R*R)*(((f[1]-ppm[v])*(f[1]+ppm[v]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]))+(f[2]-f[1])/(f[2]-f[0]));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*((f[1]-f[0])/(f[2]-f[0])+((ppm[v+1]-f[1])*(2.0*f[2]-ppm[v+1]-f[1]))/((f[2]-f[0])*(f[2]-f[1])));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[2])){
                                    bin[v]=bin[v]+1./(R*R*R);
                                }
                                else if((ppm[v]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*((ppm[v+1]-ppm[v])*(2.0*f[2]-ppm[v+1]-ppm[v]))/((f[2]-f[0])*(f[2]-f[1]));
                                }
                                else
                                    break;

                            }
                        }
                    }
                }
            }
        }
    }
    free(freq);
}

void simulate_static(int TD, vector<float>& bin, vector<float>& ppm, float diso, float span, float skew, float chi, float eta, float S, float vs, float alpha, float beta, float gamma, int fast) {
    //function to calculate a static spectrum
    int N = 32, i, j, k, v, quadrant;
    float l, m, n, ct, st, cp, sp, R, f[3];
    float spacing = ppm[1] - ppm[0];
    int Np = N + 1;
    float* freq;
    freq = (float*)malloc(4 * Np * Np * sizeof(float));

    float d22 = skew * span / 3. + diso;
    float d33 = (3. * diso - d22 - span) / 2.;
    float d11 = 3. * diso - d22 - d33;
    float ca = cos(alpha);
    float sa = sin(alpha);
    float cb = cos(beta);
    float sb = sin(beta);
    float cg = cos(gamma);
    float sg = sin(gamma);

    for (j = 0; j <= N; j++) {
        for (i = 0; (i + j) <= N; i++) {
            for (quadrant = 0; quadrant < 4; quadrant++) {
                k = N - abs(i) - abs(j);
                l = abs(i) / sqrt(i * i + j * j + k * k);
                m = abs(j) / sqrt(i * i + j * j + k * k);
                n = k / sqrt(i * i + j * j + k * k);

                ct = n;
                st = sqrt(1. - n * n);

                switch (quadrant) {
                    case 0:
                        cp = -l / st;
                        sp = m / st;
                        break;

                    case 1:
                        cp = -l / st;
                        sp = -m / st;
                        break;

                    case 2:
                        cp = l / st;
                        sp = m / st;
                        break;

                    case 3:
                        cp = l / st;
                        sp = -m / st;
                        break;
                }

                if (cp < -1.) {
                    cp = -1;
                }
                if (cp > 1.) {
                    cp = 1;
                }
                if (ct == 1.) {
                    cp = 0.;
                    st = 0.;
                }
                if (ct == 0.) {
                    st = 1.;
                }
                if (cp == 1.) {
                    sp = 0.;
                }

                float c2p = cp * cp - sp * sp;

                freq[quadrant * Np * Np + abs(i) * Np + abs(j)] = vQ_stat(diso, chi, eta, vs, S, ct, c2p) + vCSA_stat(d11, d22, d33, ct, st, cp, sp, ca, sa, cb, sb, cg, sg);
            }
        }
    }
    int a;
    for (j = 1; j <= N; j++) {
        for (i = 0; i <= (N - abs(j)); i++) {
            int jump;
            for (quadrant = 0; quadrant < 4; quadrant++) {
                for (a = 0; a < 2; a++) {
                    jump = 1;
                    if (a == 0) {
                        f[0] = freq[quadrant * Np * Np + Np * i + j];
                        f[1] = freq[quadrant * Np * Np + Np * i + j - 1];
                        f[2] = freq[quadrant * Np * Np + Np * (i + 1) + j - 1];
                    }

                    else if (a == 1 && (i + j) < (N)) {
                        f[0] = freq[quadrant * Np * Np + Np * i + abs(j)];
                        f[1] = freq[quadrant * Np * Np + Np * (i + 1) + j - 1];
                        f[2] = freq[quadrant * Np * Np + Np * (i + 1) + j];
                    } else
                        jump = 0;

                    k = N - (i) - (j);
                    R = sqrt(i * i + j * j + k * k) / N;
                    qsort(f, 3, sizeof(float), floatcomp);  // in ascending order

                    if (jump) {
                        int v1 = (int)((f[0] - ppm[0]) / spacing);
                        int v3 = (int)((f[2] - ppm[0]) / spacing);
                        int v2 = (int)((f[1] - ppm[0]) / spacing);
                        int vstart = v1;
                        int vend = v3;
                        int vmid = v2;
                        if (vstart < 0) {
                            vstart = 0;
                            v1 = -100e6;
                        }
                        if (vend > TD - 1) {
                            vend = TD - 1;
                            v3 = 100e6;
                        }
                        if (vmid < 0)
                            vmid = 0;
                        if (vmid > TD - 1)
                            vmid = TD - 1;

                        if(fast){
                            for (v = vstart; v < vmid; v++) {
                                bin[v] = bin[v] + (v1!=v2)*(v3!=v1)* (float)(v - v1) / (v2 - v1) / (R * R * R * (v3 - v1));
                            }
                            for (v = vmid; v < vend; v++) {
                                bin[v] = bin[v] + (v3!=v2)*(v3!=v1)*(float)(v3 - v) / (v3 - v2) / (R * R * R * (v3 - v1));
                            }
                        }
                        else{
                            for (v=vstart; v<TD-1; v++){
                                if((f[0]>ppm[TD-1]))
                                    break;

                                else if((f[2]<ppm[0]))
                                    break;

                                if((ppm[v]>f[0]) && (ppm[v+1]<=f[1])){
                                    bin[v]=bin[v]+1./(R*R*R)*((ppm[v+1]-ppm[v])*(ppm[v]+ppm[v+1]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[0]) && (ppm[v+1]<=f[1])){
                                    bin[v]=bin[v]+1./(R*R*R)*((ppm[v+1]-f[0])*(ppm[v+1]-f[0]))/((f[2]-f[0])*(f[1]-f[0]));
                                }
                                else if((ppm[v]>f[1]) && (ppm[v]<=f[2]) && (ppm[v+1]>f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*((f[2]-ppm[v])*(f[2]-ppm[v]))/((f[2]-f[0])*(f[2]-f[1]));
                                }
                                else if((ppm[v]>f[0]) && (ppm[v]<=f[1]) && (ppm[v+1]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*(((f[1]-ppm[v])*(f[1]+ppm[v]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]))+((ppm[v+1]-f[1])*(2.0*f[2]-ppm[v+1]-f[1]))/((f[2]-f[0])*(f[2]-f[1])));
                                }
                                else if((ppm[v]>f[0]) && (ppm[v]<=f[1]) && (ppm[v+1]>f[2])){
                                bin[v]=bin[v]+1./(R*R*R)*(((f[1]-ppm[v])*(f[1]+ppm[v]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]))+(f[2]-f[1])/(f[2]-f[0]));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*((f[1]-f[0])/(f[2]-f[0])+((ppm[v+1]-f[1])*(2.0*f[2]-ppm[v+1]-f[1]))/((f[2]-f[0])*(f[2]-f[1])));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[2])){
                                    bin[v]=bin[v]+1./(R*R*R);
                                }
                                else if((ppm[v]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+1./(R*R*R)*((ppm[v+1]-ppm[v])*(2.0*f[2]-ppm[v+1]-ppm[v]))/((f[2]-f[0])*(f[2]-f[1]));
                                }
                                else
                                    break;

                            }
                        }
                    }
                }
            }
        }
    }
    free(freq);
}

void simulate_static_sat(int TD, vector<float>& bin, vector<float>& ppm, float diso, float span, float skew, float chi, float eta, float S, float vs, float alpha, float beta, float gamma, int fast) {
    //function to calculate a static spectrum
    int N = 32, i, j, k, v, quadrant;
    float l, m, n, ct, st, cp, sp, R, f[3];
    float spacing = ppm[1] - ppm[0];
    int Np = N + 1;
    float* freq;
    freq = (float*)malloc(4 * Np * Np * sizeof(float));

    float d22 = skew * span / 3. + diso;
    float d33 = (3. * diso - d22 - span) / 2.;
    float d11 = 3. * diso - d22 - d33;
    float ca = cos(alpha);
    float sa = sin(alpha);
    float cb = cos(beta);
    float sb = sin(beta);
    float cg = cos(gamma);
    float sg = sin(gamma);

    float transition;
    for(transition = S; transition > -S; transition--){
    float moment = (S*(S+1.)-transition*(transition-1.));

    for (j = 0; j <= N; j++) {
        for (i = 0; (i + j) <= N; i++) {
            for (quadrant = 0; quadrant < 4; quadrant++) {
                k = N - abs(i) - abs(j);
                l = abs(i) / sqrt(i * i + j * j + k * k);
                m = abs(j) / sqrt(i * i + j * j + k * k);
                n = k / sqrt(i * i + j * j + k * k);

                ct = n;
                st = sqrt(1. - n * n);

                switch (quadrant) {
                    case 0:
                        cp = -l / st;
                        sp = m / st;
                        break;

                    case 1:
                        cp = -l / st;
                        sp = -m / st;
                        break;

                    case 2:
                        cp = l / st;
                        sp = m / st;
                        break;

                    case 3:
                        cp = l / st;
                        sp = -m / st;
                        break;
                }

                if (cp < -1.) {
                    cp = -1;
                }
                if (cp > 1.) {
                    cp = 1;
                }
                if (ct == 1.) {
                    cp = 0.;
                    st = 0.;
                }
                if (ct == 0.) {
                    st = 1.;
                }
                if (cp == 1.) {
                    sp = 0.;
                }

                float c2p = cp * cp - sp * sp;

                if(transition==0.5)
                    freq[quadrant * Np * Np + abs(i) * Np + abs(j)] = vQ_stat(diso, chi, eta, vs, S, ct, c2p) + vCSA_stat(d11, d22, d33, ct, st, cp, sp, ca, sa, cb, sb, cg, sg);
                else
                    freq[quadrant * Np * Np + abs(i) * Np + abs(j)] = vQ_sat(diso,chi,eta,vs,S,transition,ct,st,c2p) + vCSA_stat(d11, d22, d33, ct, st, cp, sp, ca, sa, cb, sb, cg, sg);
            }
        }
    }
    int a;
    for (j = 1; j <= N; j++) {
        for (i = 0; i <= (N - abs(j)); i++) {
            int jump;
            for (quadrant = 0; quadrant < 4; quadrant++) {
                for (a = 0; a < 2; a++) {
                    jump = 1;
                    if (a == 0) {
                        f[0] = freq[quadrant * Np * Np + Np * i + j];
                        f[1] = freq[quadrant * Np * Np + Np * i + j - 1];
                        f[2] = freq[quadrant * Np * Np + Np * (i + 1) + j - 1];
                    }

                    else if (a == 1 && (i + j) < (N)) {
                        f[0] = freq[quadrant * Np * Np + Np * i + abs(j)];
                        f[1] = freq[quadrant * Np * Np + Np * (i + 1) + j - 1];
                        f[2] = freq[quadrant * Np * Np + Np * (i + 1) + j];
                    } else
                        jump = 0;

                    k = N - (i) - (j);
                    R = sqrt(i * i + j * j + k * k) / N;
                    qsort(f, 3, sizeof(float), floatcomp);  // in ascending order

                    if (jump) {
                        int v1 = (int)((f[0] - ppm[0]) / spacing);
                        int v3 = (int)((f[2] - ppm[0]) / spacing);
                        int v2 = (int)((f[1] - ppm[0]) / spacing);
                        int vstart = v1;
                        int vend = v3;
                        int vmid = v2;
                        if (vstart < 0) {
                            vstart = 0;
                            v1 = -100e6;
                        }
                        if (vend > TD - 1) {
                            vend = TD - 1;
                            v3 = 100e6;
                        }
                        if (vmid < 0)
                            vmid = 0;
                        if (vmid > TD - 1)
                            vmid = TD - 1;

                        if(fast){
                            for (v = vstart; v < vmid; v++) {
                                bin[v] = bin[v] + moment*(v1!=v2)*(v3!=v1)* (float)(v - v1) / (v2 - v1) / (R * R * R * (v3 - v1));
                            }
                            for (v = vmid; v < vend; v++) {
                                bin[v] = bin[v] + moment*(v3!=v2)*(v3!=v1)*(float)(v3 - v) / (v3 - v2) / (R * R * R * (v3 - v1));
                            }
                        }
                        else{
                            for (v=vstart; v<TD-1; v++){
                                if((f[0]>ppm[TD-1]))
                                    break;

                                else if((f[2]<ppm[0]))
                                    break;

                                if((ppm[v]>f[0]) && (ppm[v+1]<=f[1])){
                                    bin[v]=bin[v]+moment/(R*R*R)*((ppm[v+1]-ppm[v])*(ppm[v]+ppm[v+1]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[0]) && (ppm[v+1]<=f[1])){
                                    bin[v]=bin[v]+moment/(R*R*R)*((ppm[v+1]-f[0])*(ppm[v+1]-f[0]))/((f[2]-f[0])*(f[1]-f[0]));
                                }
                                else if((ppm[v]>f[1]) && (ppm[v]<=f[2]) && (ppm[v+1]>f[2])){
                                    bin[v]=bin[v]+moment/(R*R*R)*((f[2]-ppm[v])*(f[2]-ppm[v]))/((f[2]-f[0])*(f[2]-f[1]));
                                }
                                else if((ppm[v]>f[0]) && (ppm[v]<=f[1]) && (ppm[v+1]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+moment/(R*R*R)*(((f[1]-ppm[v])*(f[1]+ppm[v]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]))+((ppm[v+1]-f[1])*(2.0*f[2]-ppm[v+1]-f[1]))/((f[2]-f[0])*(f[2]-f[1])));
                                }
                                else if((ppm[v]>f[0]) && (ppm[v]<=f[1]) && (ppm[v+1]>f[2])){
                                bin[v]=bin[v]+moment/(R*R*R)*(((f[1]-ppm[v])*(f[1]+ppm[v]-2.0*f[0]))/((f[2]-f[0])*(f[1]-f[0]))+(f[2]-f[1])/(f[2]-f[0]));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+moment/(R*R*R)*((f[1]-f[0])/(f[2]-f[0])+((ppm[v+1]-f[1])*(2.0*f[2]-ppm[v+1]-f[1]))/((f[2]-f[0])*(f[2]-f[1])));
                                }
                                else if((ppm[v]<=f[0]) && (ppm[v+1]>f[2])){
                                    bin[v]=bin[v]+moment/(R*R*R);
                                }
                                else if((ppm[v]>f[1]) && (ppm[v+1]<=f[2])){
                                    bin[v]=bin[v]+moment/(R*R*R)*((ppm[v+1]-ppm[v])*(2.0*f[2]-ppm[v+1]-ppm[v]))/((f[2]-f[0])*(f[2]-f[1]));
                                }
                                else
                                    break;

                            }
                        }
                    }
                }
            }
        }
    }

    }//transition loop
    free(freq);
}

float spectrum_RMSD(int TD, int MAS, int sat, int fast, int sites, vector<float>& bin, vector<float>& ppm, float* diso, float* span, float* skew, float* chi, float* eta, float S, float vs, float* alpha, float* beta, float* gamma,  vector<float>& LB, vector<float>& GB, float width,vector<float>& intensity) {
    //Compares a simulation to the experiment and returns the RMSD
    //The spectrum simulation and broadening functions are applied within this function
    vector<float> bin_sim(TD, 0.);
    vector<float> bin_temp(TD, 0.);
    int i, j;
    float RMSD = 0.;

    for (i = 0; i < sites; i++) {
        if (MAS) {
            simulate_MAS(TD, bin_temp, ppm, diso[i], chi[i], eta[i], S, vs,fast);
            Broadening(LB[i], GB[i], bin_temp, width, TD);
        }
        else if(sat){
            simulate_static_sat(TD, bin_temp, ppm, diso[i], span[i], skew[i], chi[i], eta[i], S, vs, alpha[i], beta[i], gamma[i],fast);
            Broadening(LB[i], GB[i], bin_temp, width, TD);
        }
        else {
            simulate_static(TD, bin_temp, ppm, diso[i], span[i], skew[i], chi[i], eta[i], S, vs, alpha[i], beta[i], gamma[i],fast);
            Broadening(LB[i], GB[i], bin_temp, width, TD);
        }
        for (j = 0; j < TD; j++) {
            bin_sim[j] = bin_sim[j] + bin_temp[j] * intensity[i];
        }
    }
    for (i = 0; i < TD; i++) {
        RMSD = RMSD + pow(bin[i] - bin_sim[i], 2.);
    }
    return 1000. * sqrt(RMSD) / TD;
}


void create_params(struct par* parameters, int nspec, vector<int>& TD, vector<int>& MAS, vector<int>& sat,int fast, int sites, vector< vector<float> >& bin, vector< vector<float> >& ppm, float S, vector<float>& vs, float* width){
    int i,j;
    parameters->nspec=nspec;
    parameters->sites=sites;
    parameters->S=S;
    parameters->ppm.resize(nspec);
    parameters->bin.resize(nspec);

    parameters->TD.resize(nspec);
    parameters->MAS.resize(nspec);
    parameters->sat.resize(nspec);
    parameters->vs.resize(nspec);
    parameters->width.resize(nspec);

    for(i=0;i<nspec;i++){
        parameters->TD[i]=TD[i];
        parameters->MAS[i]=MAS[i];
        parameters->sat[i]=sat[i];
        parameters->vs[i]=vs[i];
        parameters->width[i]=width[i];
        parameters->ppm[i].resize(TD[i],0.);
        parameters->bin[i].resize(TD[i],0.);
        for(j=0;j<TD[i];j++){
            parameters->ppm[i][j]=ppm[i][j];
            parameters->bin[i][j]=bin[i][j];
        }
    }
}

void arrange_variables(gsl_vector *var, int nspec, int sites, float* diso, float* span, float* skew, float* chi, float* eta, float* alpha, float* beta, float* gamma, vector< vector<float> >& LB, vector< vector<float> >&GB, vector< vector<float> >& intensity){
    int el=0, i, j;

    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, diso[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, span[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, skew[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, chi[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, eta[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, alpha[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, beta[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        gsl_vector_set(var, el, gamma[i]);
        el++;
    }
    for(i=0;i<sites;i++){
        for(j=0;j<nspec;j++){
            gsl_vector_set(var, el, LB[j][i]);
            el++;
        }
    }
    for(i=0;i<sites;i++){
        for(j=0;j<nspec;j++){
            gsl_vector_set(var, el, GB[j][i]);
            el++;
        }
    }
    for(i=0;i<sites;i++){
        for(j=0;j<nspec;j++){
            gsl_vector_set(var, el, intensity[j][i]);
            el++;
        }
    }
}

void disentangle_variables(const gsl_vector *var, int nspec, int sites, float* diso, float* span, float* skew, float* chi, float* eta, float* alpha, float* beta, float* gamma, vector< vector<float> >& LB, vector< vector<float> >&GB, vector< vector<float> >& intensity){
    int el=0, i, j;

    for(i=0;i<sites;i++){
        diso[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        span[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        skew[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        chi[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        eta[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        alpha[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        beta[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        gamma[i]=gsl_vector_get(var,el);
        el++;
    }
    for(i=0;i<sites;i++){
        for(j=0;j<nspec;j++){
            LB[j][i]=gsl_vector_get(var,el);
            el++;
        }
    }
    for(i=0;i<sites;i++){
        for(j=0;j<nspec;j++){
            GB[j][i]=gsl_vector_get(var,el);
            el++;
        }
    }
    for(i=0;i<sites;i++){
        for(j=0;j<nspec;j++){
            intensity[j][i]=gsl_vector_get(var,el);
            el++;
        }
    }
}

void set_step_size(gsl_vector *var, gsl_vector *ss, int sites, int nspec){
    int i,j,el=0;

    for(i=0;i<sites;i++){//diso
        gsl_vector_set(ss, el, 0.5);
        el++;
    }
    for(i=0;i<sites;i++){//span
        gsl_vector_set(ss, el, 0.05*gsl_vector_get(var,el));
        el++;
    }
    for(i=0;i<sites;i++){//skew
        gsl_vector_set(ss, el, 0.05);
        el++;
    }
    for(i=0;i<sites;i++){//chi
        gsl_vector_set(ss, el, 0.05*gsl_vector_get(var,el));
        el++;
    }
    for(i=0;i<sites;i++){//eta
        gsl_vector_set(ss, el, 0.05);
        el++;
    }
    for(i=0;i<sites;i++){//alpha
        gsl_vector_set(ss, el, 0.1);
        el++;
    }
    for(i=0;i<sites;i++){//beta
        gsl_vector_set(ss, el, 0.1);
        el++;
    }
    for(i=0;i<sites;i++){//gamma
        gsl_vector_set(ss, el, 0.1);
        el++;
    }
    for(i=0;i<sites;i++){//LB
        for(j=0;j<nspec;j++){
            gsl_vector_set(ss, el, 0.05*gsl_vector_get(var,el));
            el++;
        }
    }
    for(i=0;i<sites;i++){//GB
        for(j=0;j<nspec;j++){
            gsl_vector_set(ss, el, 0.05*gsl_vector_get(var,el));
            el++;
        }
    }
    for(i=0;i<sites;i++){//intensity
        for(j=0;j<nspec;j++){
            gsl_vector_set(ss, el, 0.03*gsl_vector_get(var,el));
            el++;
        }
    }

}

double total_RMSD(const gsl_vector* var, void* params){
    struct par *parameters = (struct par*) params;
    int k, sites=parameters->sites, nspec=parameters->nspec;
    float diso[sites], span[sites], skew[sites], chi[sites], eta[sites], alpha[sites], beta[sites], gamma[sites];
    vector<vector<float> > GB(nspec,vector<float>(sites, 0.));
    vector<vector<float> > LB(nspec,vector<float>(sites, 0.));
    vector<vector<float> > intensity(nspec,vector<float>(sites, 0.));

    disentangle_variables(var,nspec,sites,diso,span,skew,chi,eta,alpha,beta,gamma,LB,GB,intensity);

    for(k=0;k<sites;k++){
        if((fabs(skew[k])>1.)||(eta[k]>1.0)||(eta[k]<0.0)||(span[k]<0.0))
            return 1000000000000.;
    }

    float RMSD = 0.;
     for (k = 0; k < nspec; k++) {
        RMSD += spectrum_RMSD(parameters->TD[k], parameters->MAS[k], parameters->sat[k], parameters->fast, sites, parameters->bin[k], parameters->ppm[k], diso, span, skew, chi, eta, parameters->S, parameters->vs[k], alpha, beta, gamma, LB[k], GB[k], parameters->width[k], intensity[k]);
     }
    return RMSD;
}

void write_fits(float RMSD_min, int nspec, vector<int>& TD, vector<int>& MAS, vector<int>& sat, int sites,  vector< vector<float> >& bin,  vector< vector<float> >& ppm, float* diso, float* span, float* skew, float* chi, float* eta, float S, vector<float>& vs, float* alpha, float* beta, float* gamma,vector< vector<float> >& LB,vector< vector<float> >& GB,vector< vector<float> >& intensity, float* width){
    //This function creates a comma-delimited file with the experimental spectra and their simulations
    //These can be plotted in a spreadsheet program to assess the quality of the fit
    int i, j, k, TDmax=0;
    FILE *fp;
    char filename[64];
    vector< vector<float> > bin_sim(nspec,vector<float>(1, 0.));
    vector< vector<float> > bin_temp(nspec,vector<float>(1, 0.));

    for (i=0;i<nspec;i++){
        if(TD[i]>TDmax)
            TDmax=TD[i];
    }

    for(i=0;i<nspec;i++){
        bin_sim[i].resize(TD[i]);
        bin_temp[i].resize(TD[i]);

        for (j = 0; j < sites; j++) {
            if (MAS[i]) {
                simulate_MAS(TD[i], bin_temp[i], ppm[i], diso[j], chi[j], eta[j], S, vs[i],0);
                Broadening(LB[i][j], GB[i][j], bin_temp[i], width[i], TD[i]);
            }
            else if (sat[i]) {
                simulate_static_sat(TD[i], bin_temp[i], ppm[i], diso[j], span[j], skew[j], chi[j], eta[j], S, vs[i], alpha[j], beta[j], gamma[j],0);
                Broadening(LB[i][j], GB[i][j], bin_temp[i], width[i], TD[i]);
            }
            else {
                simulate_static(TD[i], bin_temp[i], ppm[i], diso[j], span[j], skew[j], chi[j], eta[j], S, vs[i], alpha[j], beta[j], gamma[j],0);
                Broadening(LB[i][j], GB[i][j], bin_temp[i], width[i], TD[i]);
            }
            for (k = 0; k < TD[i]; k++) {
                bin_sim[i][k] = bin_sim[i][k] + bin_temp[i][k] * intensity[i][j];
            }
        }
    }

    sprintf(filename,"%.5f_RMSD_fits.csv",RMSD_min);
    fp=fopen(filename,"w");

    for(j=0;j<nspec;j++){
        if(MAS[j])
            fprintf(fp,"%.1f MHz MAS,,,,",vs[j]/1.e6);
        else
            fprintf(fp,"%.1f MHz static,,,,",vs[j]/1.e6);
    }
    fprintf(fp,"\n");

    for(j=0;j<nspec;j++){
            fprintf(fp,"ppm,expt,sim,,");
    }
    fprintf(fp,"\n");

    for(i=0;i<TDmax;i++){
        for(j=0;j<nspec;j++){
            if(i<TD[j])
                fprintf(fp,"%f,%f,%f,,",ppm[j][i],bin[j][i],bin_sim[j][i]);
            else
                fprintf(fp,",,,,");
        }
        fprintf(fp,"\n");
    }

    fclose(fp);

}
