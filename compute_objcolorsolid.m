function [XYZ,V,hull]=compute_objcolorsolid(cmf,illum,wavelengths, alpha)

if nargin<4
    alpha=0.5
end


% generate samples on ocs surface using reparametrized values for a more even
% distribution of points

range=linspace(0,1,30);
N=numel(range)^2;

rOPT=zeros(N,2);

i=0;
for L1=range
    for L2=range
        i=i+1;
        rOPT(i,:)=[L1,L2];
    end
end



% convert to normal wl in nm
OPT=deparametrize(rOPT, illum, wavelengths);

% convert to XYZ

XYZ=ntr2lms(OPT,illum, cmf, wavelengths);

XYZ(XYZ<0)=0;

% draw the convex hull

% V=draw_conv_hull_color(XYZ,alpha);
[hull V]=convhulln(XYZ);
