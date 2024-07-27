function [v_OCS,OCS_EE_radii,normalizedOCS]=InertiaTensorBased_OCS_normalization(OCS1)

[OCS_id_face]=convhulln(OCS1);

%moving the OCS to the origin
[V_OCS, CR_OCS1, VRR_OCS, IRR_OCS]=polhedrn(OCS1(:,1),OCS1(:,2),OCS1(:,3),OCS_id_face);
for channel=1:3
    OCS(:,channel)=OCS1(:,channel)-CR_OCS1(channel);
end

%Principal axes of the OCS
[V_OCS, CR_OCS, VRR_OCS, IRR_OCS]=polhedrn(OCS(:,1),OCS(:,2),OCS(:,3),OCS_id_face);
[u_OCS, s_OCS, v_OCS]=svd(IRR_OCS);


%rotate the OCS
rotated_OCS=OCS*v_OCS;


%% finding the equivalent ellipsoid radii
b=rotated_OCS;
[b_id_face]=convhulln(b);
[V_b, CR_b, VRR_b, IRR_b]=polhedrn(b(:,1),b(:,2),b(:,3),b_id_face);

A=abs(IRR_b(1,1));
B=abs(IRR_b(2,2));
C=abs(IRR_b(3,3));

X2=C+B-A;
Y2=C+A-B;
Z2=B+A-C;

bCal=nthroot((15*(Y2 ^2)/(8*pi*sqrt(X2)*sqrt(Z2))),5);
aCal=nthroot((15*(X2 ^2)/(8*pi*sqrt(Y2)*sqrt(Z2))),5);
cCal=nthroot((15*(Z2 ^2)/(8*pi*sqrt(X2)*sqrt(Y2))),5);
k=sort([aCal bCal cCal]);
OCS_EE_radii=[k(1) 0 0;...
    0,k(2),0;...
    0 0 k(3)];

normalizedOCS=rotated_OCS*inv(OCS_EE_radii);