function [pS1,pS2,SS1,ind1,ind2] = phase_denoising(cS)
% phase_denoising - phase de-noising based on thresholds for spatial source phase
% Input: 
    % cS: the corrected spatial map estimate [1 x V]
% Output: 
    % pS: the corrected and denoised spatial map estimate [1 x V]
%%
pS1=zeros(size(cS));
pS2=zeros(size(cS)); SS1 = zeros(size(cS));
out_p=angle(cS);
ind1=find(out_p<=(pi/4)&out_p>=(-pi/4));
SS1(ind1) = 1;
pS1(ind1)=cS(ind1);

out_p=angle(cS);
ind2=find(out_p<(-(3*pi/4))|out_p>(3*pi/4));
pS2(ind2)=cS(ind2);
end