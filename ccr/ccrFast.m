function [D]= ccrFast(Up,Sp,Vp,Ut,St,Vt,dp,dt,dcca)
% function to make predictions using CCA-BP 
% version Fast: updates from cca.m - only generates slopes
% Written by M. Osman (osmanm@mit.edu) May 2019, modified from function
%   written by J. Smerdon (Columbia Univ.) as presented in Smerdon et al., 2010, J. Clim.
%
% Input parameters:
% Up / Sp / Vp - SVD decomposed matrices of proxy matrix P; i.e., [Up,Sp,Vp] = svd(P,'econ'); 
% Ut / St / Vt - SVD decomposed matrices of climate matrix T; i.e., [Ut,St,Vt] = svd(T,'econ'); 
% Note: P and T should have the same row dimension
% dp and dt are reduction dimensions for P and T respectively, 
% per Barnett and Preisendorfer's CCA version
% dcca is the number of CCA modes to be used for prediction
% Note: dcca cannot exceed min(dp,dt)
% 
% Output parameters:
% D - prediction operator

F = Vt(:,1:dt)'*Vp(:,1:dp);
[Uf, Sf, Vf] = svd(F);

D = Ut(:,1:dt)*St(1:dt,1:dt)*Uf(1:dt,1:dcca)*Sf(1:dcca,1:dcca)*Vf(1:dp,1:dcca)'*diag(1.0./diag(Sp(1:dp,1:dp)))*Up(:,1:dp)';
