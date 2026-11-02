function [b,cS,cA] = energymax_T1(A_time,S1,Sm2)
% energymax_T1 - Phase de-ambiguity based on maximizing the energy of the TC real part.
% Input: 
    % A_time: the time course estimate [T x 1] T is the number of time points
    % S1: the spatial map estimate [1 x V] V is the number of in-brain voxels
	% Sm2: the spatial reference network with dimension 1 x V
% Output: 
    % b: rotation angle of spatial map
    % cS: the corrected spatial map estimate [1 x V]
    % cA: the corrected time course estimate [T x 1]
%%
i=sqrt(-1); [p,q]=size(S1);
Re=real(A_time)';
Im=imag(A_time)';
A0=[Re;Im];
[V, Lambda] = icatb_v_pca(A0',1,2,0);
V=circshift(V,1);
AA=V*A0;
Ar=AA(1,:)';
Ai=AA(2,:)';
cA=Ar+i*Ai;
b=acos(V(1,1));
A_new1=A_time.*exp(i*b);
if ~eq((round(A_new1(1,1).*10000))./10000,(round(cA(1,1).*10000))./10000)
    b = -b;
    A_new1=A_time.*exp(i*b);
end
cS=S1.*exp(i*(-b));
temp=real(cS);
temp1=corrcoef(temp,Sm2);
corA=temp1(1,2);
if (corA>0)
    cA=cA;
    cS=cS;
else
    cA=-cA;
    cS=-cS;
end
