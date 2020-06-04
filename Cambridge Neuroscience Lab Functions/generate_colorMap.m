function tintColMap=generate_colorMap
%Function generates the number of colours that are used to draw the
%ratemaps

tintColMap=zeros(10, 3);

tintColMap=[1 1 1; 0 0 0.6; 0 0 0.8; 0 0.5 1; 0.2 0.7 0.9;...
    0.3 0.9 0.6;  0.4 1 0.3; 0.9 1 0.1; 1 0.8889 0; 1 0.7778 0; 1 0 0];

%colormap(tintColMap); colorbar;
end