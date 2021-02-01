function [Yprime, Gprime, skipped, deltaCalc] = TSSAT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transdelta, skipCalc, doGroup)
    deltaCalc=[];
    skipped=false;
        
    if skipCalc==false && doGroup==false
        %we are performing ITSS
        if size(transdelta,2)~=1
            error = "unexpected delta length"
            error()
        end  

        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);

        GsetminY=((G-Y)>0);
        XsetminG=((X-G)>0);
        
        if Y(xor)==1 && GsetminY(xtar)==1
            transarg = [transX transdelta];
            Yprime = FRS(G,(Sigma_c+Sigma_u)>0,transarg,Y);
            Gprime = G;
            
        elseif G(xor)==1 && XsetminG(xtar)==1 && Sigma_u(sigma)==1
            transarg = transX;
                
            %prune transarg
                i=1;
                while i<=size(transarg,2)
                    if transarg(1,i)==xor
                        transarg(:,i)=[];
                        continue;
                    end
                    i=i+1;
                end
                
            arg1 = X0;
            arg1(xor)=0;
            arg2 = Xm;
            arg2(xor)=0;
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,arg1,arg2);
            
        elseif XsetminG(xor)==1 && xor~=xtar
            transarg = [transX transdelta];
            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm)>0);
            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,arg1,arg2);
            
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        %skip FRS, BRS and SS calls
        if size(transdelta,2)~=1
            error = "unexpected delta length"
            error()
        end
        
        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
        
        GsetminY=((G-Y)>0);
        XsetminG=((X-G)>0);
        
        if Y(xor)==1 && GsetminY(xtar)==1
            skipped = true; Yprime = Y; Gprime=G;
            
        elseif G(xor)==1 && XsetminG(xtar)==1 && Sigma_u(sigma)==1
            skipped = true; Yprime = Y; Gprime=G;
            
        elseif XsetminG(xor)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
            
        else
            Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
        
        GsetminY=((G-Y)>0);
        XsetminG=((X-G)>0);
        
        Gr1=[];
        Gr2=[];
        Gr3=[];
        for trans=1:size(transdelta,2)
            td=(transdelta(:,trans));
            xor=td(1);
            sigma=td(2);
            xtar=td(3);
            
            if Y(xor)==1 && GsetminY(xtar)==1
                Gr1=[Gr1 td];
            elseif G(xor)==1 && XsetminG(xtar)==1 && Sigma_u(sigma)==1
                Gr2=[Gr2 td];
            elseif XsetminG(xor)==1 && xor~=xtar
                Gr3=[Gr3 td];
            else
                0;
            end
        end
        
        if sum(Gr1)>0
            transarg = [transX Gr1];
            Yprime = FRS(G,(Sigma_c+Sigma_u)>0,transarg,Y);
            Gprime = G;
            deltaCalc=Gr1;
            
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
                i=1;
                while i<=size(transarg,2)
                    if transarg(1,i)==xor
                        transarg(:,i)=[];
                        continue
                    end
                    i=i+1;
                end
            end
            
            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,arg1,arg2);
            deltaCalc=Gr2;
            
        elseif sum(Gr3)>0
            
            transarg = [transX Gr3];
            arg1 = ((Y+X0)>0);
            arg2 = ((Y+Xm)>0);
            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,arg1,arg2);
            deltaCalc=Gr3;
            
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end