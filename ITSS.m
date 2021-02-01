function [X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin] = ITSS(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transmin, transplus, X0plus, X0min, Xmplus, Xmmin, skipCalc, doGroup)

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
    [Y, G, skipped, ~] = TSSAT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transdelta, skipCalc, doGroup);
    if skipped==false
        transX = [transX transdelta];
    else
        transplus_skipped = [transplus_skipped transdelta];
    end
end
transplus = transplus_skipped; %this means it is empty if none were skipped, which we expect

transmin_skipped = [];
for ti = 1:size(transmin,2)
    transdelta = transmin(:,ti);
    [Y, G, skipped, ~] = TSSRT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transdelta, skipCalc, doGroup);
    if skipped==false
        edgfrmxor = find(transX(1,:)==transdelta(1));
        for j=1:length(edgfrmxor)
            t = transX(:,edgfrmxor(j));
            if t(2)==transdelta(2) && t(3)==transdelta(3)
                transX(:,edgfrmxor(j))=[];
                break; %this assumes this exact transition only exists once in transX
            end
        end
    else
        transmin_skipped = [transmin_skipped transdelta];
    end
end
transmin = transmin_skipped; %this means it is empty if none were skipped, which we expect

end