function MLsynth(foldermodelname,modelname)
% clear all; close all; clc

set_n_trans_trav(0);

% foldermodelname = 'TL_1';
% modelname = 'TL_base';
% modelname = 'TL_variant';

filename_MLsys = ['MLsys/' foldermodelname '/' modelname '_MLsys.mat'];

load(['MLsys/' foldermodelname '/' modelname '_MLsys.mat'])

[Y, G] = DLSS(X, Sigma_c, Sigma_u, transX, X0, Xm);

save(['MLsys/' foldermodelname '/' modelname '_sup_MLsys.mat'],'Y','G')

fprintf(['Traversed over ' num2str(get_n_trans_trav) ' transitions\n']);
