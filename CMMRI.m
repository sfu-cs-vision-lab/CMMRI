close all
clear

%loading camera sensor spectral sensitivities
wavelengths=[400:720]';
load('Canon500D.mat');
SS_sensors= interp1(400:10:720',Canon500D,wavelengths);
  
%loading human cone cmf
[zz,zzwl]=observer('2');
cmf_XYZ(:,1)= interp1(zzwl,zz(:,1),wavelengths);
cmf_XYZ(:,2)= interp1(zzwl,zz(:,2),wavelengths);
cmf_XYZ(:,3)= interp1(zzwl,zz(:,3),wavelengths);

%loading light
illum_info=xlsread('D65_6500k.xlsx');
illum_D65=interp1(illum_info(:,1),illum_info(:,2),wavelengths);
    
%defining the color mechanisms
maxY=illum_D65(:)'*SS_sensors(:,2);
SS_sensors=SS_sensors*100/maxY;
cm_sensor=ColorMechanism(illum_D65, SS_sensors,  wavelengths,1);

maxY=illum_D65(:)'*cmf_XYZ(:,2);
cmf_XYZ=cmf_XYZ*100/maxY;
cm_XYZ=ColorMechanism(illum_D65, cmf_XYZ,  wavelengths,1);

%Compute the OCS of the second color mechanism
[ OCS,OCS_vol,OCS_id_face]=compute_objcolorsolid(cmf_XYZ, illum_D65, wavelengths, 0.3);

%sample color to compute MMV, here we chose grey
XYZ_grey=(illum_D65'*SS_sensors)./2;
nPoints=1000; % the number of boundary points on the MMV. The higher this value, the more accurate estimation of the MMV size

%computing the MMV using Logvinenko et al. algorithm
ntrs=calculate_mmv(XYZ_grey, cm_sensor, nPoints);
MMV=cm_XYZ.ntr2lms(ntrs);
[~,MMV_vol]=convhulln(MMV);

%normalize the MMV based on OCS
[v_OCS1,OCS_EE_radii1,normalizedOCS1]=InertiaTensorBased_OCS_normalization(OCS);

%compute the MMV dimension using its equivalent ellipsoid
[MMV_Dim1,IRR1,normalizedMMV1]=InertiaTensorBased_MMV_normalization(MMV,v_OCS1,OCS_EE_radii1);

%display the result
disp('The average of MMV dimensions')
mean(MMV_Dim1)

