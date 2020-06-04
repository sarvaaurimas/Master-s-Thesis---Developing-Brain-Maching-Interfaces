function Draw_RMap(rmap)
%draws the specified rate map 'rmap'
pixels_per_m1 = 350;
bin_size1=ceil(pixels_per_m1*0.025);

NrColors=10;
tintColMap=generate_colorMap;



%normalized xy coordinates, i.e. counted from 1 to the length inpixels of
%the environment

tmp_value=min(min(rmap))-(max(max(rmap))-min(min(rmap))+1)/NrColors;
tmp_rmap=reshape(rmap,1,size(rmap,1)*size(rmap,2));
tmp_rmap(isnan(tmp_rmap))=tmp_value;
rmap2=reshape(tmp_rmap, size(rmap,1),size(rmap,2));



set(gcf, 'color', [1 1 1]);
imagesc(rmap2);
colormap(tintColMap);
set(gca,'YDir','reverse');
axis equal;
xticklabels([0, 20, 40, 60, 80, 100, 120]*2.5);
yticklabels([1:11]*10*2.5);
%axis off;
%axis tight;
xlabel('x, cm','FontSize',14);
ylabel('y, cm','FontSize',14);

end