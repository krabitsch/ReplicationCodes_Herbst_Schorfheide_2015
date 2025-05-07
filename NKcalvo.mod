// Dynare implementation of the Basic New Keynesian Example Model used in Herbst and Schorfheide (2015): Bayesian Estimation of DSGE Models
// NK model with Calvo pricing
// Katrin Rabitsch, 2018

// Endogenous variables
var Y g R PI dA S K F N C Z obs_dY obs_PI obs_R;

// Innovations
varexo epsZ epsG epsR;

// PARAMETERS
parameters tau epsil phi nu theta rho_R psi1 psi2 xiH gammaQ rA piA Zbar rho_z sigma_z Gfrac rho_g sigma_g PItarget
           C_ss N_ss S_ss R_ss K_ss F_ss PI_ss dA_ss G_ss Y_ss g_ss; // PItarget gama betta

// fixed parameters:
xiH        = 1;
phi        = 0;             // inverse Frisch elast. of labor (=0 implies disutility linear in labor)
epsil      = 6;             // elasticity between varieties, 6 -> 20 percent (net) markup, 11 -> 10 percent markup (not separately identified from theta and phi) 
nu         = 1/epsil;       // 1/epsil -> removes monopolistic competition distortion
Zbar       = 1;             // stst value of tech
Gfrac      = 0.2;           // fraction of output that falls on gov expend

// parameters to be estimated:
tau        = 2;             // coeff. of RRA
theta      = 0.66;          // Calvo parameter: 0.75 -> 4 quarters sticky

psi1       = 2.50;          // Taylor rule coefficient on inflation
psi2       = 0.25;          // Taylor rule coefficient on output gap
rA         = 0.30;          // annualized real interest rate
piA        = 4.00;          // annualized inflation
gammaQ     = 0.5;           // quarterly growth rate

rho_R      = 0.8;           // interest rate smoothing in Taylor
rho_g      = 0.8;           // persistence of gov shock
rho_z      = 0.9;           // persistence of tech shock
sigma_R    = 0.004;         // vol. of mon. pol. shock
sigma_g    = 0.01;          // vol. of gov. shock
sigma_z    = 0.005;         // vol. of tech shock

PItarget   = 1;             // Schorfheide's log-linearized model is around 0-stst-inflation (even if piA>1)
//PItarget = piA   /400+1;  // Inflation target (and stst inflation)
//gama     = gammaQ/100+1;  // gross growth rate
//betta    = 1/(1+rA/400);  // discount factor
//kappa = (((1-theta)*(1-theta*betta))/theta) * ((tau+phi)/1)



model;

//#PItarget  = piA   /400+1;  // Inflation target (and stst inflation)
#gama      = gammaQ/100+1;  // gross growth rate
#betta     = 1/(1+rA/400);  // discount factor
#Rsteady = (gama/betta)*PItarget;
#Yflex = ( (1-nu)*(epsil/(epsil-1)) * (1/(xiH*(1-Gfrac)^tau)) )^(1/(tau+phi));

// Euler equation:
exp(C)^(-tau) = (betta*exp(C(+1))^(-tau))*(exp(R)/exp(PI(+1))) * (1/exp(dA(+1)));

// optimal price setting, Calvo: 
(((1-theta*exp(PI)^(epsil-1))/(1-theta))^(1/(1-epsil))) = exp(K)/exp(F);

// auxiliary equations for Calvo price setting:
exp(K) =  exp(Y) * (1 - nu)*(epsil/(epsil-1))*exp(C)^(tau)*xiH*(exp(N)^phi) + betta*theta* (exp(C(+1))/exp(C))^(-tau) * (exp(PI(+1))^ epsil   ) * exp(K(+1));
exp(F) =  exp(Y)                                                            + betta*theta* (exp(C(+1))/exp(C))^(-tau) * (exp(PI(+1))^(epsil-1)) * exp(F(+1));

// price dispersion: 
exp(S) = ( (1-theta)*(((1-theta*(exp(PI))^(epsil-1))/(1-theta))^(epsil/(epsil-1))) + theta*(exp(PI)^epsil)*exp(S(-1)) );

// resource constraint: 
exp(C)*exp(g) = (exp(Y)/exp(S));

// production function:
exp(Y) = exp(N);

// Taylor rule (in the Schorfheide version (where the output component is in terms of Y-G ): 
//exp(R)/Rsteady = (exp(R(-1))/Rsteady)^rho_R * ( (exp(PI)/PItarget)^psi1  * (exp(Y)/Y_ss)^psi2 )^(1-rho_R) * exp(epsR);
exp(R)/Rsteady = (exp(R(-1))/Rsteady)^rho_R * ( (exp(PI)/PItarget)^psi1  * ( (exp(Y)/Y_ss) / (exp(g)/g_ss) )^psi2 )^(1-rho_R) * exp(epsR);

// economy's growth rate
exp(dA) = gama*exp(Z);

// Exogenous processes: 
(Z)  = rho_z  * (Z(-1))                          + epsZ;
(g)  = rho_g  * (g(-1)) + (1-rho_g) * log(g_ss)  + epsG;

// Observables: 
//obs_dY = gammaQ + Y - Y(-1) + Z;
//obs_PI = piA + 4*PI;
//obs_R  = piA + rA + 4*gammaQ + 4*R;
obs_dY = gammaQ + log(exp(Y)) - log(exp(Y(-1))) + log(exp(Z)/Zbar) ;
obs_PI = piA + 4*log(exp(PI)/PI_ss);
obs_R  = piA + rA + 4*gammaQ + 4*log(exp(R)/R_ss);

end;



resid(1);
steady;
check;

shocks;
var epsZ;
stderr sigma_z;
var epsG;
stderr sigma_g;
var epsR;
stderr sigma_R;
end;

stoch_simul(order=1) Y C g PI R S Z dA N; 


//stoch_simul(order=2,periods=5000,nograph,pruning); 

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
theta,        0.66   , 0.0001  , 0.999 , BETA_PDF,     0.66,    0.2;
psi1,         2.50   , 1.0001  , 10    , GAMMA_PDF,     1.5,    0.25;
psi2,         0.25   , 0.0001  , 5     , GAMMA_PDF,     0.5,    0.25;
rA,           0.30   , 0.01    , 2     , GAMMA_PDF,     0.3,    0.2;
piA,             4   , 0.001   , 12    , GAMMA_PDF,       4,       2;
gammaQ,       0.5    , 0.001   , 2     , NORMAL_PDF,    0.5,     0.2;
rho_R,        0.8    ,   0.0   , 0.999 , BETA_PDF,      0.8,     0.15;
rho_g,        0.8    ,   0.001 , 0.999 , BETA_PDF,      0.8,     0.15;
rho_z,        0.9    ,   0.001 , 0.999 , BETA_PDF,      0.8,     0.15;
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

estimation(optim=('MaxIter',200),datafile=NKmodel_Schorfheide_data,mode_compute=1,mode_check,mh_replic=50000,mh_nblocks=2,mh_jscale=0.50);

