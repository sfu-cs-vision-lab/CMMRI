function [ OPT ] = deparametrize( rOPT, cmf , wavelengths)
% DEPARAMETRIZE takes Logvinenko's reparametrized wavelengths and returns
% regular wavelengths. rOPT can be of any dimension.

% create the deparametrization function
sigma=cumtrapz(sqrt(sum(cmf.^2,2)));
sigma=sigma/max(sigma);    
invomegaspline=spline(sigma, wavelengths);
invomega=@(wl) ppval(invomegaspline,wl);

OPT=rOPT*0;

OPT(:)=invomega(rOPT(:));