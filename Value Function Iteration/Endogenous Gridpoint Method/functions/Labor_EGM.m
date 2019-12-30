% Labor_EGM.m is a routine for computing labor satisfying eq. (17) in the
% neoclassical growth model with elastic labor supply, considered in the 
% article "Envelope Condition Method versus Endogenous Grid Method for Solving 
% Dynamic Programming Problems" by Lilia Maliar and Serguei Maliar,  Economics 
% Letters (2013), 120, 262-266 (henceforth, MM, 2013).
% -------------------------------------------------------------------------
% Inputs:    "dvf1" is next-period derivative of value functions;
%            "k0" is current capital;
%            "z0" is current productivity;
%            "A" is the normalizing constant in production;
%            "alpha" is the share of capital in production;
%            "gam, nu, B" are the utility-function parameters;
%            "delta" is the depreciation rate;
%            "beta" is the discount factor; 
%            "x" is unknown labor choice
%
% Output:    "dif_FOC" is the difference between the right and left hand 
%            sides of the FOC 
% -------------------------------------------------------------------------
% Copyright © 2011-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [dif_FOC] = Labor_EGM(x,nu,gam,alpha,delta,beta,B,A,k1,z0,dvf1)

dif_FOC=k1+(beta*dvf1)^(-1/gam)-(1-delta)*(B*(1-x)^-nu/z0/A/(1-alpha)/dvf1/beta)^(1/alpha)*x-B*(1-x)^-nu/(1-alpha)/dvf1/beta*x;
                            % Evaluate the difference between the left and right
                            % sides of the FOC 
