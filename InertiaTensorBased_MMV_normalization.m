function [MMV_Dim,IRR_b,normalizedMMV]=InertiaTensorBased_MMV_normalization(MMV,v_OCS,OCS_EE_radii)


%normalizing the MMV
normalizedMMV=MMV*v_OCS*inv(OCS_EE_radii);

%% Compute the normalized MMV dimension
[MMV_id_face]=convhulln(normalizedMMV);

%move MMV to the origin
[V_MMV, CR_MMV, VRR_MMV, IRR_MMV]=polhedrn(normalizedMMV(:,1),normalizedMMV(:,2),normalizedMMV(:,3),MMV_id_face);
for channel=1:3
    New_MMV(:,channel)=normalizedMMV(:,channel)-CR_MMV(channel);
end

%Principal axes of translated MMV
[V_MMV, CR_MMV, VRR_MMV, IRR_MMV]=polhedrn(New_MMV(:,1),New_MMV(:,2),New_MMV(:,3),MMV_id_face);
[u_MMV, s_MMV, v_MMV]=svd(IRR_MMV);

%rotate the MMV
rotated_MMV=New_MMV*v_MMV;


%inertia tensor of the normalized MMV
[b_id_face]=convhulln(rotated_MMV);
[V_b, CR_b, VRR_b, IRR_b]=polhedrn(rotated_MMV(:,1),rotated_MMV(:,2),rotated_MMV(:,3),b_id_face);


%% 
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
MMV_Dim=k;

