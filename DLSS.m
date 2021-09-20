%% DLSS (Ouedrago synth)
function [Y, G] = DLSS(X, Sigma_c, Sigma_u, transK, X0, Xm)
    Xk1=X;
    while(true)
        Nk = BRS(Xk1,(Sigma_c+Sigma_u)>0,transK,Xk1.*Xm);
        Bk = BRS(Xk1,Sigma_u,transK,(X.*(Xk1-Nk))>0);
        Xk2 = X.*((Xk1-Bk))>0;
        if isequal(Xk2,Xk1)
            G = Xk2;
            break;
        else
            Xk1 = Xk2;
        end
    end
    Y = FRS(G,Sigma_c+Sigma_u,transK,X0.*G);
end