function [X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin] = ITSS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transminindices, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup)

if doGroup==true
    error();
end

for xd=1:length(X0plus)
    if X0plus(xd)==1
        [Y, G, skipped, ~] = TSSAIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            X0(xd)=1;
            X0plus(xd)=0;
        end
    end
end
for xd=1:length(X0min)
    if X0min(xd)==1
        [Y, G, skipped, ~] = TSSRIS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            X0(xd)=0;
            X0min(xd)=0;
        end
    end
end

for xd=1:length(Xmplus)
    if Xmplus(xd)==1
        [Y, G, skipped, ~] = TSSAMS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            Xm(xd)=1;
            Xmplus(xd)=0;
        end
    end
end
for xd=1:length(Xmmin)
    if Xmmin(xd)==1
        [Y, G, skipped, ~] = TSSRMS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, xd, skipCalc, doGroup);
        if skipped==false
            Xm(xd)=0;
            Xmmin(xd)=0;
        end
    end
end

transplus_skipped = [];
for ti = 1:size(transplus,2)
    transdelta = transplus(:,ti);
    [Y, G, skipped, ~, ~] = TSSAT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transdelta, skipCalc, doGroup);
    if skipped==false
        transX = [transX transdelta];
    else
        transplus_skipped = [transplus_skipped transdelta];
    end
end
transplus = transplus_skipped;

transmin_skipped = [];
for ti = transminindices
    [Y, G, skipped, ~] = TSSRT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, ti, skipCalc, doGroup);
    if skipped==false
        transX(:,ti)=[0;0;0];
    else
        transmin_skipped = [transmin_skipped ti];
    end
end
transminindices = transmin_skipped;

end