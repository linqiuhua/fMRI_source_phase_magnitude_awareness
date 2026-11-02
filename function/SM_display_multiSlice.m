function SM_display_multiSlice(SS,thr,comid)
cd SMshow
colortables
anatomy = icatb_spm_read_vols(icatb_spm_vol('rsingle_subj_T1.img'));
out=zeros(53,63,46);
out(:)=SS(:);
figure_FontSize=24;

if comid==1;
    max1 = roundn(maxN(out(:,:,[18:22,25,27,28,29])),-2);
    min1 = thr;
   for threshold=thr
           max(reshape(out(1:1:end,1:1:end,[18:22,25,27,28,29]),[],1));
           out=make_composite2(out(1:1:end,1:1:end,[18:22,25,27,28,29]),anatomy(1:1:end,1:1:end,[18:22,25,27,28,29]),threshold,'nearest'); %%ÏÔÊ¾1·ù
           if min1>=0,
               axis image;colormap(hot);colorbar('YTick',[0  255.5],'YTickLabel',[min1 max1'],'FontSize',figure_FontSize,'FontName','Arial Unicode MS')
           else
               axis image;colormap(coldhot);colorbar('YTick',[-0.5  255.5],'YTickLabel',[min1 max1'],'FontSize',figure_FontSize,'FontName','Arial Unicode MS')
           end   
    end
end
if comid==2;
    max1 = roundn(maxN(out(:,:,[14:22])),-2);
    min1 = thr;

   for threshold=thr
       disp(roundn(max(max(max(out(:,:,[14:22])))),-2));
       out=make_composite2(out(1:1:end,1:1:end,[14:22]),anatomy(1:1:end,1:1:end,[14:22]),threshold,'nearest'); %%ÏÔÊ¾1·ù
       if min1>=0,
           axis image;colormap(hot);colorbar('YTick',[-0.5  255.5],'YTickLabel',[min1 max1'],'FontSize',figure_FontSize,'FontName','Arial Unicode MS')
       else
           axis image;colormap(coldhot);colorbar('YTick',[-0.5  255.5],'YTickLabel',[min1 max1'],'FontSize',figure_FontSize,'FontName','Arial Unicode MS')
       end   
   end
end

if comid==3;
    max1 = roundn(maxN(out(:,:,[18:26])),-2);
     min1 = thr;
   for threshold=thr
         disp(roundn(max(max(max(out(:,:,[18:26])))),-2));
         out=make_composite2(out(1:1:end,1:1:end,[18:26]),anatomy(1:1:end,1:1:end,[18:26]),threshold,'nearest'); %%ÏÔÊ¾1·ù
         if min1>=0
             axis image;colormap(hot);colorbar('YTick',[-0.5  255.5],'YTickLabel',[min1 max1'],'FontSize',figure_FontSize,'FontName','Arial Unicode MS')
        else
             axis image;colormap(coldhot);colorbar('YTick',[-0.5  255.5],'YTickLabel',[min1 max1'],'FontSize',figure_FontSize,'FontName','Arial Unicode MS')
        end
    end
end

cd ..