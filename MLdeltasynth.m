function time = MLdeltasynth(foldermodelname,basemodelname,variantmodelname,doITSS)
filename_P_allinit_statespace = ['CIFreport/' foldermodelname '/' basemodelname '_allinit_statespace.cif'];
filename_MLsys = ['MLsys/' foldermodelname '/' basemodelname '_MLsys.mat'];
filename_MLsup = ['MLsys/' foldermodelname '/' basemodelname '_sup_MLsys.mat'];
filename_modeldelta = ['MLsys/' foldermodelname '/' basemodelname '_to_' variantmodelname '_modeldelta.mat'];

load(filename_MLsys)
load(filename_MLsup)
load(filename_modeldelta)

    X = ((X+Xplus)>0); Sigma_c = ((Sigma_c+Sigma_cplus)>0); Sigma_u = ((Sigma_u+Sigma_uplus)>0);
    Xplus = zeros(1,length(X)); Sigma_cplus = zeros(1,length(Sigma_c)); Sigma_uplus = zeros(1,length(Sigma_u)); 

    tic

%% Only perform these lines for ITSS
if doITSS==true
    skipCalc=false;
    doGroup=false;
    
    [X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin] = ITSS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup);
     
%% These lines for GTSS
else
       
    while true
        skipCalc=true;
        doGroup=false;
        [X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin] = ITSS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup);
        DeltaSize = size(transplus,2)+size(transminindices,2)+sum(X0plus)+sum(X0min)+sum(Xmplus)+sum(Xmmin);
        if DeltaSize==0
            break
        end

        skipCalc=false;
        doGroup=true;
        [X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin] = GTSS(X, Sigma_c, Sigma_u, transX, adj, adjrev, adj_urev, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup);
        DeltaSize = size(transplus,2)+size(transminindices,2)+sum(X0plus)+sum(X0min)+sum(Xmplus)+sum(Xmmin);
        if DeltaSize==0
            break
        end
        
    end
    
end

    time = toc;

    X=((X-Xmin)>0); Sigma_c=((Sigma_c-Sigma_cmin)>0); Sigma_u=((Sigma_u-Sigma_umin)>0);
    Xmin = zeros(1,length(X)); Sigma_cmin = zeros(1,length(Sigma_c)); Sigma_umin = zeros(1,length(Sigma_u)); 

save(['MLsys/' foldermodelname '/' variantmodelname '_deltasup_MLsys.mat'], 'X', 'Sigma_c', 'Sigma_u', 'transX', 'X0', 'Xm', 'Y', 'G')

end
