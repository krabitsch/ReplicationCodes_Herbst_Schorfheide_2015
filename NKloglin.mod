// NK model of Schorfheide slides, Bayesian estimation

// Endogenous variables
var Y PI R G Z C obs_dY obs_PI obs_R;

// Innovations
varexo epsZ epsG epsR;

// PARAMETERS
parameters tau kappa psi1 psi2 rA piA gammaQ rho_R rho_g rho_z sigma_R sigma_g sigma_z betta;

tau     = 2.00;
//kappa   = 0.175406324044845;
kappa   = 0.350812648089690; //0.2;
//kappa = 3.030302990318019e+08;
psi1    = 2.50;
psi2    = 0.25;
rA      = 0.30;
piA     = 4.00;
gammaQ  = 0.5;
rho_R   = 0.8;
rho_g   = 0.8;
rho_z   = 0.9;
sigma_R = 0.4;
sigma_g = 1.00;
sigma_z = 0.5;

betta = 1/(1+rA/400);


model(linear);
#gama      = gammaQ/100+1;  // gross growth rate

// Euler equation:
Y = Y(+1) - (1/tau) * ( R - PI(+1) - Z(+1) ) + G - G(+1);

// Phillips-curve: 
//PI = betta * PI(+1) + kappa * ( Y - G );
PI = betta * PI(+1) + kappa * ( Y - G );

// Taylor rule: 
R = rho_R * R(-1) + (1-rho_R) * psi1 * PI + (1-rho_R) * psi2 * ( Y - G ) + epsR; 
//R = rho_R * R(-1) + (1-rho_R) * psi1 * PI + (1-rho_R) * psi2 * ( Y  ) + epsR; 

// Resource constraint:
Y  = C + G;

// Exogenous processes: 
G  = rho_g * G(-1) + epsG;
Z  = rho_z * Z(-1) + epsZ;

// Observables: 
obs_dY = gammaQ + Y - Y(-1) + Z;
obs_PI = piA + 4*PI;
obs_R  = piA + rA + 4*gammaQ + 4*R;

end;



resid(1);
steady;
check;

shocks;
var epsZ; stderr sigma_z;
var epsG; stderr sigma_g;
var epsR; stderr sigma_R;
end;

stoch_simul(order=2) Y C G PI R Z; 

//stoch_simul(order=2,periods=500,nograph); 



estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
//  1) tau is GAMMA with mean 2 and st.d. 0.50
//  2) kappa is UNIFORM with [0,1]
//  3) psi1 is GAMMA with mean 1.50 and st.d 0.25
//  4) psi2 is GAMMA with mean 0.50 and st.d.0.25
//  5) rA is GAMMA with mean 0.50 and st.d. 0.50
//  6) piA is GAMMA with mean 7.00 and variance 2.00
//  7) gammaQ is Normal with mean 0.40 and variance 0.20
//  8) rho_R is UNIFORM with [0,1)
//  9) rho_g is UNIFORM with [0,1)
// 10) rho_z is UNIFORM with [0,1)
// 11) sigma_R is InvGAMMA with mean 0.40 and st.d. 4.00
// 12) sigma_g is InvGAMMA with mean 1.00 and st.d. 4.00
// 13) sigma_z is InvGAMMA with mean 0.50 and st.d. 4.00


tau,          2.00   , 0.5     , 5     , GAMMA_PDF,       2,    0.5;
kappa,        0.20   , 0.0001  , 0.9   , BETA_PDF,      0.2,    0.1;
psi1,         2.50   , 1.0001  , 10    , GAMMA_PDF,     1.5,    0.25;
psi2,         0.25   , 0.0001  , 5     , GAMMA_PDF,     0.5,    0.25;
rA,           0.30   , 0.01    , 2     , GAMMA_PDF,     0.3,    0.2;
piA,             4   , 0.001   , 12    , GAMMA_PDF,       4,       2;
gammaQ,       0.5    , 0.001   , 2     , NORMAL_PDF,    0.5,     0.2;
rho_R,        0.8    ,   0.0   , 0.999 , BETA_PDF,      0.8,     0.15;
rho_g,        0.8    ,   0.001 , 0.999 , BETA_PDF,      0.8,     0.15;
rho_z,        0.9    ,   0.001 , 0.999 , BETA_PDF,      0.9,     0.10;
stderr epsR,  0.40   ,  0.00001, 10    , INV_GAMMA_PDF, 0.4,     4;
stderr epsG,  1.00   ,  0.00001, 10    , INV_GAMMA_PDF, 1.0,     4;
stderr epsZ,  0.5    ,  0.00001, 10    , INV_GAMMA_PDF, 0.5,     4;


/*
tau,          2.00   , 0.5     , 5     ;
kappa,        0.20   , 0.0001  , 0.9   ;
psi1,         2.50   , 1.0001  , 10    ;
psi2,         0.25   , 0.0001  , 5     ;
rA,           0.30   , 0.01    , 2     ;
piA,             4   , 0.001   , 12    ;
gammaQ,       0.5    , 0.001   , 2     ;
rho_R,        0.8    ,   0.0   , 0.999 ;
rho_g,        0.8    ,   0.001 , 0.999 ;
rho_z,        0.9    ,   0.001 , 0.999 ;
stderr epsR,  0.40   ,  0.00001, 10    ;
stderr epsG,  1.00   ,  0.00001, 10    ;
stderr epsZ,  0.5    ,  0.00001, 10    ;
*/

end;

varobs obs_dY obs_PI obs_R;

estimation(optim=('MaxIter',200),datafile=NKmodel_Schorfheide_data,mode_compute=1,mode_check,mh_replic=10000,mh_nblocks=2,mh_jscale=0.50);


