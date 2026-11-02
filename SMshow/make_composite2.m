function out = make_composite(func,anat,thresh,ker,gridsize)
%func and anat must have the same number of slices
%function out = make_composite(func,anat,thresh,ker)

func = squeeze(func);%移除单唯独

if (~exist('ker','var')),%如何调整图片
   ker = 'bilinear';
end;
   
if (~exist('gridsize','var')),
   gridsize = -999;
end;

colortables

xdim = size(anat,1);ydim = size(anat,2);zdim = size(anat,3);
xdimf = size(func,1);ydimf = size(func,2);zdimf = size(func,3);

%resize functional map to be larger (same as the template)
clear func2;
for j = 1:zdim,
   func2(:,:,j) = imresize(func(:,:,j),[xdim ydim],ker);%调整图像 维数和加权插值
end;
func = func2;
%%
%scale functional to range from 129-256
if (minN(func)>=0),%N维最小  恒为false 去else
%%
   %positive only
   v = [];
   ind = find(func>thresh);if (length(ind)<2), v = 0; end; 
   func = (128+(func.*(func>thresh)-minN([func(func>thresh);v]))* ...
      127/maxN(func-minN([func(func>thresh);v]))).*(func>thresh);
   clear gray;
   h = hot(384);
%    aa = [gray(256);hot];
    aa = [gray(256);h([120:end],:)];   %%%%may chang color-lin
   
%    load hot.mat
%    aa = [gray(256);hot];   %%%%may chang color-dan
   %aa = [gray(256);imresize(cold(60:end,:),[256 3])];
   
   %gr = gray(1000);
   %aa = [gray(256);gr(end:-1:(end-255),:)];
else,
   %positive and negative
   [Y I]=max(flatrow(abs(func)));%取出func中绝对值的最大值和索引
   Y = func(I);%取出最大值

   func(1)=abs(Y);%+max
   func(end)=-abs(Y);%-max
   
   func = (128+(func.*(abs(func)>thresh)-minN(func))* ...
      127/maxN(func-minN(func))).*(abs(func)>thresh);%归一化到128-256 详见本子 2016.1.8
   func(1) = 0;func(end)=0;%加入的用于计算的值移除
    
    
   clear gray;
%    aa = [gray(256);hot]; %RGB的值
     aa = [gray(256);coldhot]; %RGB的值  
%    func = (128+(func.*(func>thresh)-minN(func))* ...
%       127/maxN(func-minN(func))).*(func>thresh);
%    clear gray;
%    aa = [gray(256);coldhot];
   
end;

%scale anatomic to range from 1-128
anat = (anat-minN(anat))*127/maxN(anat-minN(anat));%归一化到1-128

%combine them
out = anat.*(func==0)+func;%把我们的和脑结构的结合


%pick which slices to display...
%out=out(:,:,:);
% %%%%%%%%%%%%%%%添加小图--------修改
% outall=anat.*(func==0)+func;%%--------改动
% out=zeros(53,63,9);
% out(:,:,1)=outall(:,:,4);
% out(:,:,2)=outall(:,:,8);
% out(:,:,3)=outall(:,:,9);
% 
% out(:,:,4)=outall(:,:,16);
% out(:,:,5)=outall(:,:,18);
% out(:,:,6)=outall(:,:,19);
% 
% out(:,:,7)=outall(:,:,20);
% out(:,:,8)=outall(:,:,21);
% out(:,:,9)=outall(:,:,22);

%%%%%%%%%%%%%%%%%%%%%%%%%
a=create_2d(out,size(out,1),size(out,2),size(out,3),1,size(out,3),1,1,gridsize,0,256,0);%生成信号颜色如图

% colormap(aa(1:2:end,:))
subimage(round(a')+1,aa(1:2:end,:));axis off;%加入脑结构颜色图
%set(gcf,'position',[22   506   518   429]);
%print(gcf,'-dtiff',['nocrave_' num2str(jj,'%.2d')]);


