function [norm_pos_map, norm_rate_map] = normalize_maps( pos_map, rate_map, pos_sampl_rate)    
% Normalises the dwellmaps so it is seconds spent per bin
% Normalises the ratemaps so it is spikes per second per bin

norm_pos_map = pos_map ./ pos_sampl_rate;
norm_rate_map = rate_map ./ norm_pos_map;
