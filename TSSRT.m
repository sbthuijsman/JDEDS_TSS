function [Yprime, Gprime, skipped, deltaCalc] = TSSRT(X, Sigma_c, Sigma_u, transX, X0, Xm, Y, G, transdelta, skipCalc, doGroup)
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
        
        if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
            transarg = transX;
            
            edgfrmxor = find(transarg(1,:)==xor);

            for j=1:size(edgfrmxor,2)
                t = transarg(:,edgfrmxor(j));
                if t(2)==sigma && t(3)==xtar
                    transarg(:,edgfrmxor(j))=[];
                    break; %this assumes this exact transition only exists once in transX
                end
            end
            
            if size(transarg,2)~=size(transX,2)-1
                error = "unexpected transarg length"
                error();
            end        

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,X0,Xm);
            
        elseif GsetminY(xor)==1 && G(xtar)==1 && xor~=xtar
            transarg = transX;
            
            edgfrmxor = find(transarg(1,:)==xor);

            for j=1:size(edgfrmxor,2)
                t = transarg(:,edgfrmxor(j));
                if t(2)==sigma && t(3)==xtar
                    transarg(:,edgfrmxor(j))=[];
                    break; %this assumes this exact transition only exists once in transX
                end
            end
            
            if size(transarg,2)~=size(transX,2)-1
                error = "unexpected transarg length"
                error();
            end        

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(Y+Xm)>0);
            
        elseif XsetminG(xor)==1 && XsetminG(xtar)==1 && Sigma_u(sigma)==1 && xor~=xtar
            transarg = transX;
            
            edgfrmxor = find(transarg(1,:)==xor);

            for j=1:size(edgfrmxor,2)
                t = transarg(:,edgfrmxor(j));
                if t(2)==sigma && t(3)==xtar
                    transarg(:,edgfrmxor(j))=[];
                    break; %this assumes this exact transition only exists once in transX
                end
            end
            
            if size(transarg,2)~=size(transX,2)-1
                error = "unexpected transarg length"
                error();
            end        

            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(G+Xm)>0);
            
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
        
        if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
             
        elseif GsetminY(xor)==1 && G(xtar)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;

        elseif XsetminG(xor)==1 && XsetminG(xtar)==1 && Sigma_u(sigma)==1 && xor~=xtar
            skipped = true; Yprime = Y; Gprime=G;
            
        else
            Yprime = Y; Gprime=G;
        end
        
    elseif skipCalc==false && doGroup==true
        %perform grouped FRS, BRS, or SS call
        
        xor=transdelta(1);
        sigma=transdelta(2);
        xtar=transdelta(3);

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
            
            if Y(xor)==1 && Y(xtar)==1 && xor~=xtar
                Gr1=[Gr1 td];
            elseif GsetminY(xor)==1 && G(xtar)==1 && xor~=xtar
                Gr2=[Gr2 td];
            elseif XsetminG(xor)==1 && XsetminG(xtar)==1 && Sigma_u(sigma)==1 && xor~=xtar
                Gr3=[Gr3 td];
            else
                0;
            end
        end
                
        if sum(Gr1)>0
            transarg = transX;
            
            %for all trans in group
            for ti = 1:size(Gr1,2)
                xor = Gr1(1,ti);
                sigma = Gr1(2,ti);
                xtar = Gr1(3,ti);
            
                edgfrmxor = find(transarg(1,:)==xor);

                for j=1:length(edgfrmxor)
                    t = transarg(:,edgfrmxor(j));
                    if t(2)==sigma && t(3)==xtar
                        transarg(:,edgfrmxor(j))=[];
                        break; %this assumes this exact transition only exists once in transX
                    end
                end
            end
            
%             if size(transarg,2)~=size(transX,2)-size(transarg,2)
%                 error = "unexpected transarg length"
%                 error();
%             end        

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,X0,Xm);
            deltaCalc=Gr1;
            
        elseif sum(Gr2)>0
            transarg = transX;
            
            %for all trans in group
            for ti = 1:size(Gr2,2)
                xor = Gr2(1,ti);
                sigma = Gr2(2,ti);
                xtar = Gr2(3,ti);
            
                edgfrmxor = find(transarg(1,:)==xor);

                for j=1:length(edgfrmxor)
                    t = transarg(:,edgfrmxor(j));
                    if t(2)==sigma && t(3)==xtar
                        transarg(:,edgfrmxor(j))=[];
                        break; %this assumes this exact transition only exists once in transX
                    end
                end
            end
            
%             if length(transarg)~=length(transX)-length(transarg)
%                 error = "unexpected transarg length"
%                 error();
%             end        

            [Yprime,Gprime] = DLSS(G,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(Y+Xm)>0);
            deltaCalc=Gr2;
            
        elseif sum(Gr3)>0
            transarg = transX;
            
            %for all trans in group
            for ti = 1:size(Gr3,2)
                xor = Gr3(1,ti);
                sigma = Gr3(2,ti);
                xtar = Gr3(3,ti);
            
                edgfrmxor = find(transarg(1,:)==xor);

                for j=1:length(edgfrmxor)
                    t = transarg(:,edgfrmxor(j));
                    if t(2)==sigma && t(3)==xtar
                        transarg(:,edgfrmxor(j))=[];
                        break; %this assumes this exact transition only exists once in transX
                    end
                end
            end
            
%             if length(transarg)~=length(transX)-length(transarg)
%                 error = "unexpected transarg length"
%                 error();
%             end        

            [Yprime,Gprime] = DLSS(X,Sigma_c,Sigma_u,transarg,(Y+X0)>0,(G+Xm)>0);
            deltaCalc=Gr3;
            
        else
            skipped = true; Yprime = Y; Gprime=G;
        end
        
    end
end