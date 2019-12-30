% Labor_ECM.m is a routine for computing labor satisfying eq. (18) in the
% neoclassical growth model with elastic labor supply, considered in the 
% article "Envelope Condition Method versus Endogenous Grid Method for Solving 
% Dynamic Programming Problems" by Lilia Maliar and Serguei Maliar,  Economics 
% Letters (2013), 120, 262-266 (henceforth, MM, 2013).
% -------------------------------------------------------------------------
% Inputs:    "dvf0" is the current derivative of value functions;
%            "k0" is current capital;
%            "z0" is current productivity;
%            "A" is the normalizing constant in production;
%            "alpha" is the share of capital in production;
%            "nu, B" are the utility-function parameters;
%            "delta" is the depreciation rate;
%            "x" is unknown labor choice
%
% Output:    "dif_FOC" is the difference between the right and left hand 
%            sides of the FOC 
% -------------------------------------------------------------------------
% Copyright © 2011-2016 by Lilia Maliar and Serguei Maliar. All rights 
% reserved. The code may be used, modified and redistributed under the 
% terms provided in the file "License_Agreement.txt".
% -------------------------------------------------------------------------

function [dif_FOC] = Labor_ECM(x,nu,alpha,delta,B,A,k0,z0,dvf0)

dif_FOC=dvf0*x^-alpha-B*(1-x)^-nu/z0/A/(1-alpha)/k0^alpha*(1-delta+alpha*A*z0*k0^(alpha-1)*x^(1-alpha));
                            % Evaluate the difference between the left and 
                            % right sides of the FOC 
