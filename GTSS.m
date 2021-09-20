function [X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin] = GTSS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup)

if doGroup==false
    error();
end

if sum(X0plus)>0
    [Y, G, skipped, deltaCalc] = TSSAIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, X0plus, skipCalc, doGroup);
    if skipped==false
        X0 = (X0+deltaCalc)>0;
        X0plus = (X0plus-deltaCalc)>0;
        return
    end
end

if sum(X0min)>0
    [Y, G, skipped, deltaCalc] = TSSRIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, X0min, skipCalc, doGroup);
    if skipped==false
        X0 = (X0-deltaCalc)>0;
        X0min = (X0min-deltaCalc)>0;
        return
    end
end

if sum(Xmplus)>0
    [Y, G, skipped, deltaCalc] = TSSAMS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, Xmplus, skipCalc, doGroup);
    if skipped==false
        Xm = (Xm+deltaCalc)>0;
        Xmplus = (Xmplus-deltaCalc)>0;
        return
    end
end


if sum(Xmmin)>0
    [Y, G, skipped, deltaCalc] = TSSRMS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, Xmmin, skipCalc, doGroup);
    if skipped==false
        Xm = (Xm-deltaCalc)>0;
        Xmmin = (Xmmin-deltaCalc)>0;
        return
    end
end

if size(transplus,2)>0
    [Y, G, skipped, deltaCalc, deltaCalcInd] = TSSAT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transplus, skipCalc, doGroup);
    if skipped==false
        transX = [transX deltaCalc];
        %remove deltaCalc from transplus, has been processed
        transplus(:,deltaCalcInd)=[];

        
        return
    end
end

if size(transminindices,2)>0
    [Y, G, skipped, deltaCalc] = TSSRT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transminindices, skipCalc, doGroup);
    if skipped==false
        %remove deltaCalc from transX, has been processed
        transX(:,deltaCalc)=zeros(3,size(deltaCalc,2));

        %remove deltaCalc from transmin, has been processed
        transminindices = setdiff(transminindices,deltaCalc);

        return
    end
end

end