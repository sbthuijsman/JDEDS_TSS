function [Yprime, Gprime, skipped, deltaCalc] = TSSRT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transminind, skipCalc, doGroup)
    deltaCalc=[];
    skipped=false;
    
    transdelta = transX(:,transminind);
    
    if skipCalc==false && doGroup==false
        %we are performing ITSS
        
        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
        
        if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
            
            transarg=transX;
            transarg(:,transminind)=[0;0;0];  

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,X0,Xm);
            
        elseif (G(xor)-Y(xor))==1 && G(xtar)==1 && xor~=xtar
            
            transarg=transX;
            transarg(:,transminind)=[0;0;0];

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(Y+Xm)>0);
            
        elseif (X(xor)-G(xor))==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1 && xor~=xtar
            
            transarg=transX;
            transarg(:,transminind)=[0;0;0];      

            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(G+Xm)>0);
            
        else
            Yprime = Y; Gprime=G;
        end

    elseif skipCalc==true && doGroup==false
        %skip FRS, BRS and SS calls


        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);
        
        if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
             
        elseif (G(xor)-Y(xor))==1 && G(xtar)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;

        elseif (X(xor)-G(xor))==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
            
        else
            Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
        
        Gr1=[];
        Gr2=[];
        Gr3=[];
        
        for ti = transminind
            td=transX(:,ti);
            xor=td(1);
            sigma=td(2);
            xtar=td(3);
            
            if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
                Gr1=[Gr1 ti];
            elseif (G(xor)-Y(xor))==1 && G(xtar)==1 && xor~=xtar
                Gr2=[Gr2 ti];
            elseif (X(xor)-G(xor))==1 && (X(xtar)-G(xtar))==1 && Sigma_u(sigma)==1 && xor~=xtar
                Gr3=[Gr3 ti];
            else
                0;
            end
        end
        
        if size(Gr1,2)>0
            
            transarg=transX;
            transarg(:,Gr1)=zeros(3,size(Gr1,2));   

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,X0,Xm);
            deltaCalc=Gr1;
            
        elseif size(Gr2,2)>0

            transarg=transX;
            transarg(:,Gr2)=zeros(3,size(Gr2,2));   

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(Y+Xm)>0);
            deltaCalc=Gr2;
            
        elseif size(Gr3,2)>0

            transarg=transX;
            transarg(:,Gr3)=zeros(3,size(Gr3,2));
        
            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(G+Xm)>0);
            deltaCalc=Gr3;
            
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end