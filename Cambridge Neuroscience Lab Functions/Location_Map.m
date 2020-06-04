function rmap=Location_Map(xy, bin_size, varargin)
%This function generates the position rate map in arbitrary units
%(number of visits to the bin), which you then have to devide by 
%position sampling rate to get the occupation time in sec for each bin.
%The size of the position rate map is Nr_bins x Nr_bins

%bin_size must be given in cm
%additional input represents a feature of the map that the User must specify, which defines
%whether the map will be drawn to fit tightly the trajectories (in this
%case varargin='tight', or it will be drawn having a preset size. 
%e.g. rmap=Location_Map(xy, bin_size, 230, 140) would mean that I would
%like to project my path on 230 cm (x) x 140 cm (y) plane.




Nr_pixels=bin_size; 

if length(varargin)<=1
  Nr_bins_x=ceil(max(xy(:, 1))/Nr_pixels);
  Nr_bins_y=ceil(max(xy(:, 2))/Nr_pixels);
else
  
    Nr_bins_x=ceil(varargin{1}/Nr_pixels);
    Nr_bins_y=ceil(varargin{2}/Nr_pixels);
    
end


%visitted positions expressed in bins
xy_bin(:,1)=ceil(xy(:, 1)/Nr_pixels);
xy_bin(:,2)=ceil(xy(:, 2)/Nr_pixels);

rmap=NaN(Nr_bins_y,Nr_bins_x); % Empty map 

for i=1:Nr_bins_x % for each bin in x direction
    for j=1:Nr_bins_y % for each bin in y direction
      tmpx=xy_bin(:, 1)==i; % xy is a n x 2 vector with first column x, second y
      tmpy=xy_bin(:, 2)==j;
      sum_tmp=sum(tmpx.*tmpy); % pick out only those values at that bin
      if sum_tmp~=0 %if not equal 0
        rmap(j, i)=sum_tmp;
      end
    end
end

clear xy_bin;

end