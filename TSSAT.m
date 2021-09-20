function [Yprime, Gprime, skipped, deltaCalc, deltaCalcInd] = TSSAT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transdelta, skipCalc, doGroup)
    deltaCalc=[];
    deltaCalcInd=-1;
    skipped=false;
        
    if skipCalc==false && doGroup==false
        %we are performing ITSS

        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
        
        if Y(xor)==1 && (G(xtar)-Y(xtar))==1
            transarg = [transX transdelta];
            Yprime = FRS(G,(Sigma_c+Sigma_u)>0,transarg,Y);
            Gprime = G;
            
        elseif G(xor)==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1
            transarg = transX;
                
            %prune transarg
            fnd = find(transX(1,:)==xor);
            transarg(:,fnd)=zeros(3,size(fnd,2));
                
            arg1 = X0;
            arg1(xor)=0;
            arg2 = Xm;
            arg2(xor)=0;
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,arg1,arg2);
            
        elseif (X(xor)-G(xor))==1 && xor~=xtar
            transarg = [transX transdelta];
            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm)>0);
            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,arg1,arg2);
            
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        %skip FRS, BRS and SS calls
        
        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
                
        if Y(xor)==1 && (G(xtar)-Y(xtar))==1
            skipped = true; Yprime = Y; Gprime=G;
            
        elseif G(xor)==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1
            skipped = true; Yprime = Y; Gprime=G;
            
        elseif (X(xor)-G(xor))==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
            
        else
            Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
                
        Gr1=[];
        Gr2=[];
        Gr3=[];
        Gr1ind=[];
        Gr2ind=[];
        Gr3ind=[];
        for trans=1:size(transdelta,2)
            td=(transdelta(:,trans));
            xor=td(1);
            sigma=td(2);
            xtar=td(3);
            
            if Y(xor)==1 && (G(xtar)-Y(xtar))==1
                Gr1=[Gr1 td];
                Gr1ind=[Gr1ind trans];
            elseif G(xor)==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1
                Gr2=[Gr2 td];
                Gr2ind=[Gr2ind trans];
            elseif (X(xor)-G(xor))==1 && xor~=xtar
                Gr3=[Gr3 td];
                Gr3ind=[Gr3ind trans];
            else
                0;
            end
        end
        
        if sum(Gr1)>0
            transarg = [transX Gr1];
            Yprime = FRS(G,(Sigma_c+Sigma_u)>0,transarg,Y);
            Gprime = G;
            deltaCalc=Gr1;
            deltaCalcInd=Gr1ind;
            
        elseif sum(Gr2)>0
            transarg = transX;

            arg1 = X0;
            arg2 = Xm;
            
            %for all trans in group
            for ti = 1:size(Gr2,2)
                xor = Gr2(1,ti);
                
                arg1(xor)=0;
                arg2(xor)=0;
            %prune transarg
            fnd = find(transX(1,:)==xor);
            transarg(:,fnd)=zeros(3,size(fnd,2));
            
            end
            
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,arg1,arg2);
            deltaCalc=Gr2;
            deltaCalcInd=Gr2ind;
            
        elseif sum(Gr3)>0
            
            transarg = [transX Gr3];
            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm)>0);
            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,arg1,arg2);
            deltaCalc=Gr3;
            deltaCalcInd=Gr3ind;
            
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end