function time = MLsynth(foldermodelname,modelname)

load(['MLsys/' foldermodelname '/' modelname '_MLsys.mat'])

tic

[Y, G] = DLSS(X, adj, adjrev, adj_urev, X0, Xm);

time = toc;

save(['MLsys/' foldermodelname '/' modelname '_sup_MLsys.mat'],'Y','G')

end

