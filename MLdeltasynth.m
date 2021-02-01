function MLdeltasynth(foldermodelname,basemodelname,variantmodelname,doITSS)
% clear all; close all; clc

set_n_trans_trav(0);

% foldermodelname = 'TL_1';
% basemodelname = 'TL_base';
% variantmodelname = 'TL_variant';

% doITSS = false; %GTSS if false

filename_MLsys = ['MLsys/' foldermodelname '/' basemodelname '_MLsys.mat'];
filename_MLsup = ['MLsys/' foldermodelname '/' basemodelname '_sup_MLsys.mat'];
filename_modeldelta = ['MLsys/' foldermodelname '/' basemodelname '_to_' variantmodelname '_modeldelta.mat'];

load(filename_MLsys)
load(filename_MLsup)
load(filename_modeldelta)

Delta = [sum(Xplus) sum(Xmin) sum(Sigma_cplus+Sigma_uplus) sum(Sigma_cmin+Sigma_umin) size(transplus,2) size(transmin,2) sum(X0plus) sum(X0min) sum(Xmplus) sum(Xmmin)]

%% Only perform these lines for ITSS
if doITSS==true
    skipCalc=false;
    doGroup=false;
    
    X = ((X+Xplus)>0); Sigma_c = ((Sigma_c+Sigma_cplus)>0); Sigma_u = ((Sigma_u+Sigma_uplus)>0);
    Xplus = zeros(1,length(X)); Sigma_cplus = zeros(1,length(Sigma_c)); Sigma_uplus = zeros(1,length(Sigma_u)); 
    
    [X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin] = ITSS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup);
    
    X=((X-Xmin)>0); Sigma_c=((Sigma_c-Sigma_cmin)>0); Sigma_u=((Sigma_u-Sigma_umin)>0);
    Xmin = zeros(1,length(X)); Sigma_cmin = zeros(1,length(Sigma_c)); Sigma_umin = zeros(1,length(Sigma_u)); 

    %sanity check, delta should be empty
    Delta = [sum(Xplus) sum(Xmin) sum(Sigma_cplus+Sigma_uplus) sum(Sigma_cmin+Sigma_umin) size(transplus,2) size(transmin,2) sum(X0plus) sum(X0min) sum(Xmplus) sum(Xmmin)]
    if sum(Delta)>0
        error = "unexpected deltas left"
        error();
    end

    fprintf(['Traversed over ' num2str(get_n_trans_trav) ' transitions\n']);
 
%% These lines for GTSS
else
    
    X = ((X+Xplus)>0); Sigma_c = ((Sigma_c+Sigma_cplus)>0); Sigma_u = ((Sigma_u+Sigma_uplus)>0);
    Xplus = zeros(1,length(X)); Sigma_cplus = zeros(1,length(Sigma_c)); Sigma_uplus = zeros(1,length(Sigma_u)); 
    
    while true
        skipCalc=true;
        doGroup=false;
        [X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin] = ITSS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup);
        DeltaCalc = [size(transplus,2) size(transmin,2) sum(X0plus) sum(X0min) sum(Xmplus) sum(Xmmin)]
        transitions = get_n_trans_trav
        if sum(DeltaCalc)==0
            break
        end
        skipCalc=false;
        doGroup=true;
        [X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin] = GTSS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup);
        DeltaCalc = [size(transplus,2) size(transmin,2) sum(X0plus) sum(X0min) sum(Xmplus) sum(Xmmin)]
        transitions = get_n_trans_trav
        if sum(DeltaCalc)==0
            break
        end
        
    end
    
    X=((X-Xmin)>0); Sigma_c=((Sigma_c-Sigma_cmin)>0); Sigma_u=((Sigma_u-Sigma_umin)>0);
    Xmin = zeros(1,length(X)); Sigma_cmin = zeros(1,length(Sigma_c)); Sigma_umin = zeros(1,length(Sigma_u)); 

    fprintf(['Traversed over ' num2str(get_n_trans_trav) ' transitions\n']);

end

save(['MLsys/' foldermodelname '/' variantmodelname '_deltasup_MLsys.mat'], 'X', 'Sigma_c', 'Sigma_u', 'transX', 'X0', 'Xm', 'Y', 'G')

