function time = MLsynth(foldermodelname,modelname)

load(['MLsys/' foldermodelname '/' modelname '_MLsys.mat'])

tic

[Y, G] = DLSS(X, Sigma_c, Sigma_u, transX, X0, Xm);

time = toc;

save(['MLsys/' foldermodelname '/' modelname '_sup_MLsys.mat'],'Y','G')

end

