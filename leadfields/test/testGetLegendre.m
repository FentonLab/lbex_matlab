clear all;

library = 'E:\hiren\lbex';
addpath( genpath( library ) );
% 
%   Program to test getLegendre.m against MATLAB's legendre function. 
%   While the approx. eeg leadfield algorithm does not require more than
%   3/4 terms, 
%   Test for high degree, since for direct summation of EEG leadfield
%   series, we might need to retain many terms. 
%
i=200; degree = 60;
x=[rand(i,1); -rand(i,1)];
% Our function
disp('Computing using getLegendre:');
tic, [p0,p1]=getLegendre( x, degree ); toc, ours = [p0(:,end)'; p1(:,end)'];
% Matlab's function
disp('Computing using Matlabs legendre:');
tic, tmp = legendre( degree, x, 'norm' ); toc, theirs = tmp(1:2, :);
% Errors
p0MaxError = max( abs(ours(1,:) - theirs(1,:)) )
p1MaxError = max( abs(ours(2,:) - theirs(2,:)) )

