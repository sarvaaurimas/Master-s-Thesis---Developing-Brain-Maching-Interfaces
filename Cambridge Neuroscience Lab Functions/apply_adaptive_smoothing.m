function  [smooth_pos_map, smooth_r_map]=apply_adaptive_smoothing(pos_map, spk_map, pos_sample_rate)
%Applies the adaptive smoothing on the rate map
%pos_map - unsmoothened position map
%spk_map - unsmoothened spike map

alpha=5000;

[l w]=size(pos_map);

tmp_pos_map=NaN(l+40, w+40);
tmp_spk_map=tmp_pos_map; %creating the extended matrices to place the position map and the spike map, so that the circle mask can be applied for the edges of the real map

tmp_pos_map(20:20+l-1, 20:20+w-1)=pos_map;
tmp_spk_map(20:20+l-1, 20:20+w-1)=spk_map;
tmp_smooth_r_map=NaN(size(tmp_spk_map));
tmp_smooth_pos_map=NaN(size(tmp_pos_map));

for j=20:(20+l-1)
    for i=20:(20+w-1)
      if isnan(tmp_pos_map(j, i))~=1 % if the cell has been visited
        r=1;
        ratio=10; %the initial arbitrary chosen value
        while (r<ratio)&&(j>r+1)&&(i>r+1)&&(i+r+1<=size(tmp_spk_map,2))&&(j+r+1<=size(tmp_spk_map,1)) % this makes sure the radius of filter doesn't take beyond edges
          r=r+1;
          [x, y]=meshgrid(i-r:i+r, j-r:j+r);
          c_mask=((x-i).^2+(y-j).^2)<=r^2; %a small circular mask
          map_mask=NaN(size(tmp_spk_map));
          map_mask(j-r:j+r, i-r:i+r)=c_mask;
          
          tmp1=map_mask.*tmp_pos_map;
          tmp2=map_mask.*tmp_spk_map;
          tmp1=reshape(tmp1, 1, size(tmp1, 1)*size(tmp1, 2)); %flatten the arrays
          tmp2=reshape(tmp2, 1, size(tmp2, 1)*size(tmp2, 2));
          tmp1=tmp1(isnan(tmp1)~=1); %all not NaN position bins within the circle
          tmp2=tmp2(isnan(tmp2)~=1); %all not NaN spike bins within the circle
          n=sum(tmp1); % new N value for pos map bin
          s=sum(tmp2); % new S value for spike map bin
          nr_bins = size(tmp1, 2);
          ratio=alpha/(n*sqrt(s));
          clear x y c_mask map_mask tmp1 tmp2;
        end
        tmp_smooth_r_map(j, i)=s/n*pos_sample_rate; % normalised values 
        tmp_smooth_pos_map(j, i)=n/pos_sample_rate/nr_bins;
        clear n s nr_bins;
      end
    end

end

smooth_r_map=tmp_smooth_r_map(20:20+l-1, 20:20+w-1);
smooth_pos_map=tmp_smooth_pos_map(20:20+l-1, 20:20+w-1);


end