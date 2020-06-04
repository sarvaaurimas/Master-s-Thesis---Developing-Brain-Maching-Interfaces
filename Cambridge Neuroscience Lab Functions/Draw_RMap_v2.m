function Draw_RMap_v2(rmap)
%draws the specified rate map 'rmap'

NrColors=10;
tintColMap=generate_colorMap;



%normalized xy coordinates, i.e. counted from 1 to the length inpixels of
%the environment

tmp_value=min(min(rmap))-(max(max(rmap))-min(min(rmap))+1)/NrColors;
tmp_rmap=reshape(rmap,1,size(rmap,1)*size(rmap,2));
tmp_rmap(isnan(tmp_rmap))=tmp_value;
rmap2=reshape(tmp_rmap, size(rmap,1),size(rmap,2));


set(gca, 'color', [1 1 1]);
imagesc(rmap2);
colormap(gca,tintColMap);
set(gca,'YDir','reverse');
axis equal;
%axis off;
%axis tight;


end