function [ys,params,check] = NKcalvo_steadystate(ys,exo,M_,options_)
% function [ys,params,check] = NK_baseline_steadystate(ys,exo,M_,options_)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs:
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%   - M_        [structure] Dynare model structure
%   - options   [structure] Dynare options structure
%
% Output:
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - params    [vector] vector of parameter values
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impose restrictions on parameters)

% Copyright (C) 2013-2020 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
    paramname = M_.param_names{ii};
    eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;


%%
% % for older Dynare versions (before versions 5) use this code block instead
% % of the one above

% function [ys,check] = NKcalvo_steadystate(ys,exe);
% global M_ lgy_
% 
% %% DO NOT CHANGE THIS PART.
% %%
% %% Here we load the values of the deep parameters in a loop.
% %%
% if isfield(M_,'param_nbr') == 1
%     NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
% for i = 1:NumberOfParameters                                  % Loop...
%   paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
%   eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
% end                                                           % End of the loop.  
% check = 0;
% end
% %%
%% END OF THE FIRST MODEL INDEPENDENT BLOCK.


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%

% PItarget = piA   /400+1;   
PItarget = 1;   
gama     = gammaQ/100+1;  
betta    =  1/(1+rA/400); 
Rsteady  = (gama/betta)*PItarget;
Yflex    = ( (1-nu)*(epsil/(epsil-1)) * (1/(xiH*(1-Gfrac)^tau)) )^(1/(tau+phi));


Z_ss     = Zbar;
dA_ss    = gama*Z_ss;
PI_ss    = PItarget;
R_ss     = (gama/betta)*PI_ss;

S_ss     = (1/(1-theta*(PI_ss^epsil))) *  ( (1-theta)*(((1-theta*(PI_ss)^(epsil-1))/(1-theta))^(epsil/(epsil-1))) );

% (((1-theta*PI_ss^(epsil-1))/(1-theta))^(1/(1-epsil)))  = (  K_ss/F_ss )
% K = ( Y * (1 - nu)*(epsil/(epsil-1))*C^(tau)*xiH*(N^phi) + betta*theta* (C/C)^(-tau) * (PI^ epsil   ) * K );
% F = ( Y                                                  + betta*theta* (C/C)^(-tau) * (PI^(epsil-1)) * F );
% K     = (1/(1-betta*theta* (PI_ss^ epsil   ))) * ( Y_ss * (1 - nu)*(epsil/(epsil-1))*C_ss^(tau)*xiH*(N_ss^phi)  )
% F     = (1/(1-betta*theta* (PI_ss^(epsil-1)))) * ( Y_ss                                                         )
coeff1   = (((1-theta*PI_ss^(epsil-1))/(1-theta))^(1/(1-epsil)));
% % coeff2   = (1/(1-betta*theta* dA_ss * (PI_ss^ epsil   )));
% % coeff3   = (1/(1-betta*theta* dA_ss * (PI_ss^(epsil-1))));
coeff2   = (1/(1-betta*theta* (PI_ss^ epsil   )));
coeff3   = (1/(1-betta*theta* (PI_ss^(epsil-1))));
% coeff1 = (coeff2 * (1 - nu)*(epsil/(epsil-1))*C_ss^(tau)*xiH*(N_ss^phi)) / coeff3
% coeff1*(coeff3/coeff2) = (1 - nu)*(epsil/(epsil-1)) * (C_ss)^(tau)*xiH*(N_ss^phi)
% coeff1*(coeff3/coeff2) = (1 - nu)*(epsil/(epsil-1)) * ((1-Gfrac)*Y_ss/S_ss)^(tau)*xiH*(Y_ss^phi)
% (coeff1*(coeff3/coeff2)) / ((1 - nu)*(epsil/(epsil-1)) * ((1-Gfrac)/S_ss)^(tau)*xiH) = Y_ss^(tau+phi)
Y_ss     = ( (coeff1*(coeff3/coeff2)) / ((1 - nu)*(epsil/(epsil-1)) * ((1-Gfrac)/S_ss)^(tau)*xiH) )^(1/(tau+phi));
C_ss     = (1-Gfrac)*Y_ss/S_ss;
G_ss     = Gfrac*Y_ss/S_ss;
N_ss     = Y_ss;
K_ss     = (1/(1-betta*theta* (PI_ss^ epsil   ))) * ( Y_ss * (1 - nu)*(epsil/(epsil-1))*C_ss^(tau)*xiH*(N_ss^phi)  );
F_ss     = (1/(1-betta*theta* (PI_ss^(epsil-1)))) * ( Y_ss                                                         );
Ynet_ss  = Y_ss/S_ss;
Yflex_ss = Yflex;
g_ss     = (1/(1-Gfrac));


ifprint=0;
if ifprint == 1
fprintf(' C     = %8.8f\n',C_ss );
fprintf(' N     = %8.8f\n',N_ss);
fprintf(' S     = %8.8f\n',S_ss );
fprintf(' R     = %8.8f\n',R_ss);
fprintf(' K     = %8.8f\n',K_ss );
fprintf(' F     = %8.8f\n',F_ss);
fprintf(' PI    = %8.8f\n',PI_ss );
fprintf(' G     = %8.8f\n',G_ss );
fprintf(' Y     = %8.8f\n',Y_ss);
fprintf(' Ynet  = %8.8f\n',Ynet_ss);
fprintf(' Yflex = %8.8f\n',Yflex_ss);
end




% assigning stst values to variables carried in the system:
C     = log(C_ss);
N     = log(N_ss);
S     = log(S_ss);
R     = log(R_ss);
K     = log(K_ss);
F     = log(F_ss);
PI    = log(PI_ss);
dA    = log(dA_ss);
Z     = log(Z_ss);
Y     = log(Y_ss);
G     = log(G_ss);
g     = log(g_ss);
obs_dY = gammaQ ;
obs_PI = piA ;
obs_R  = piA + rA + 4*gammaQ ;




%%
%% END OF THE MODEL SPECIFIC BLOCK.


%% DO NOT CHANGE THIS PART.
%%
params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
    eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
    varname = M_.endo_names{ii};
    eval(['ys(' int2str(ii) ') = ' varname ';']);
end


%%
% % % for older Dynare versions (before versions 5) use this code block instead
% % % of the one above
% 
% % % Here we define the steady state values of the endogenous variables of
% % % the model.
% %
% 
% for iter = 1:length(M_.params)
%   eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
% end
% 
% if isfield(M_,'param_nbr') == 1
% 
% if isfield(M_,'orig_endo_nbr') == 1
% NumberOfEndogenousVariables = M_.orig_endo_nbr;
% else
% NumberOfEndogenousVariables = M_.endo_nbr;                   % Number of endogenous variables.
% end
% ys = zeros(NumberOfEndogenousVariables,1);                   % Initialization of ys (steady state).
% for i = 1:NumberOfEndogenousVariables                        % Loop...
%   varname = deblank(M_.endo_names(i,:));                     % Get the name of endogenous variable i
%   eval(['ys(' int2str(i) ') = ' varname ';']);               % Get the steady state value of this variable.
% end                                                          % End of the loop.
% else
% ys=zeros(length(lgy_),1);
% for i = 1:length(lgy_)
%     ys(i) = eval(lgy_(i,:));
% end
% check = 0;
% end
end

