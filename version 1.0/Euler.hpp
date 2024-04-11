float Pi=3.1415926;

void Fix_Angles(float &alpha, float &beta, float &gamma){
    //Sets the Euler angles such that alpha and beta are between 0 and 90 and gamma is between 0 and 180
    //Equivalent angle sets are taken from Svenningsson and  Mueller, SSNMR, 123 (2023) 101849.
    if(beta<0.){
        beta+=Pi;
        alpha*=-1.;
    }
    if(alpha<0.)
        alpha+=Pi;
    if(gamma<0.)
        gamma+=Pi;

    float alphas[16], betas[16], gammas[16];
    int i;

    alphas[0]=alpha;         betas[0]=beta;      gammas[0]=gamma;
    alphas[1]=alpha+Pi;      betas[1]=Pi-beta;   gammas[1]=2.*Pi-gamma;
    alphas[2]=alpha+Pi;      betas[2]=Pi-beta;   gammas[2]=Pi-gamma;
    alphas[3]=alpha;         betas[3]=beta;      gammas[3]=gamma+Pi;
    alphas[4]=2.*Pi-alpha;   betas[4]=Pi-beta;   gammas[4]=gamma+Pi;
    alphas[5]=Pi-alpha;      betas[5]=beta;      gammas[5]=Pi-gamma;
    alphas[6]=Pi-alpha;      betas[6]=beta;      gammas[6]=2.*Pi-gamma;
    alphas[7]=2.*Pi-alpha;   betas[7]=Pi-beta;   gammas[7]=gamma;
    alphas[8]=Pi-alpha;      betas[8]=Pi-beta;   gammas[8]=gamma+Pi;
    alphas[9]=2.*Pi-alpha;   betas[9]=beta;      gammas[9]=Pi-gamma;
    alphas[10]=2.*Pi-alpha;  betas[10]=beta;     gammas[10]=2.*Pi-gamma;
    alphas[11]=Pi-alpha;     betas[11]=Pi-beta;  gammas[11]=gamma;
    alphas[12]=alpha+Pi;     betas[12]=beta;     gammas[12]=gamma;
    alphas[13]=alpha;        betas[13]=Pi-beta;  gammas[13]=2.*Pi-gamma;
    alphas[14]=alpha;        betas[14]=Pi-beta;  gammas[14]=Pi-gamma;
    alphas[15]=alpha+Pi;     betas[15]=beta;     gammas[15]=gamma+Pi;

    //setting angles between 0 and 2Pi
    for(i=0;i<16;i++){
        if(alphas[i]>=2.*Pi)
            alphas[i]-=2.*Pi;
        else if(alphas[i]<0.)
            alphas[i]+=2.*Pi;

        if(betas[i]>=2.*Pi)
            betas[i]-=2.*Pi;
        else if(betas[i]<0.)
            betas[i]+=2.*Pi;

        if(gammas[i]>=2.*Pi)
            gammas[i]-=2.*Pi;
        else if(gammas[i]<0.)
            gammas[i]+=2.*Pi;
    }

    //Finding the wanted equivalent angle set
    for(i=0;i<16;i++){
        if((alphas[i]<Pi/2.)&&(betas[i]<=Pi/2.)&&(gammas[i]<Pi)){
            alpha=alphas[i];  beta=betas[i];  gamma=gammas[i];
            break;
        }
    }
}

void Calculate_DMFIT_angles(float alpha, float beta, float gamma, float skew, float &phi, float &chi, float &psi){
    //Function to calculate the DMFIT/SIMPSON Euler angles using the WSOLIDS/QUEST convention.

    if(skew>=0.){
        phi = -1.*gamma;
        chi = -1.*beta;
        psi = Pi/2.-1.*alpha;
    }

    else{
        float R12= sin(alpha)*cos(beta)*cos(gamma) + cos(alpha)*sin(gamma);
        float R13=  -sin(beta)*cos(gamma);
        float R33= cos(beta);

        chi=-acos(R13);
        psi=-asin(-1.*R12/sin(chi))+Pi/2.;
        phi=-acos(R33/sin(chi));
    }
    Fix_Angles(phi,chi,psi);
}

void Calculate_WSolids_angles(float &alpha, float &beta, float &gamma, float skew, float phi, float chi, float psi){

     if(skew>=0.){
        alpha = Pi/2.-psi;
        beta=-1.*chi;
        gamma = -1.*phi;
    }

    else{
        float R12= sin(Pi/2.-psi)*cos(-chi)*cos(-phi) + cos(Pi/2.-psi)*sin(-phi);
        float R13=  -sin(-chi)*cos(-phi);
        float R23=sin(chi)*sin(phi);
        beta=-acos(R13);
        gamma=-asin(R23/sin(beta));
        alpha=asin(-1.*R12/sin(beta));

        if(psi>Pi/2.)
            gamma=asin(R23/sin(beta));

    }
    Fix_Angles(alpha,beta,gamma);
}
