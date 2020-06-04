function sacSmooth=get_smooth_SAC(spk_rmap)
%get smoothened spatial autocorrelation function(based on CB function)


  %vars=default_vars;
  gausSigma=2;
  [x,y]=meshgrid(-2:2);
  gausKern=1*(exp(-(x.^2 + y.^2)/(2*gausSigma^2)));
  gausKern=gausKern./sum(gausKern(:)); %Normalise
  clear x y

  sac = spatialAutoCorr(spk_rmap);
  
  sacSmooth=sac;
  sacSmooth(isnan(sac))=0; %Set nans to 0 to prevent nan proliferation
  
  sacSmooth=filter2(gausKern, sacSmooth);%Smooth first

end