function corellogram(x_pred, x_real)
hist3([x_pred, x_real],'Nbins', [15, 15], 'CdataMode','auto')
xlabel('X_{predicted}')
ylabel('X_{real}')
colorbar
view(2)