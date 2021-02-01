%% DLSS (Ouedrago synth)
function [Y, G] = DLSS(X, Sigma_c, Sigma_u, transK, X0, Xm)
    Xk1=X;
    iteration=0;
%      transK = prune_trans(transK,Xk1);
    while(true)
        iteration=iteration+1
%         current_n_states = sum(Xk1)
        Nk = BRS(Xk1,(Sigma_c+Sigma_u)>0,transK,Xk1.*Xm);
        Bk = BRS(Xk1,Sigma_u,transK,(X.*(Xk1-Nk))>0);
        Xk2 = X.*((Xk1-Bk))>0;
%          transK = prune_trans(transK,Xk2);
        if isequal(Xk2,Xk1)
            G = Xk2;
            break;
        else
            Xk1 = Xk2;
        end
    end
    Y = FRS(G,Sigma_c+Sigma_u,transK,X0.*G);
%     transY = prune_trans(transK,Y);
end