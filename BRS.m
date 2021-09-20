function Xplus = BRS(X,Sigma,trans,Xinit)
%reverse trans for BRS
trans = [trans(3,:)
         trans(2,:)
         trans(1,:)];

%this is same as FRS
Xplus = FRS(X,Sigma,trans,Xinit);

end