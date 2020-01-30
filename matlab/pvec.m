function c=pvec(a,b)
% PVEC Producto vectorial de dos vectores c=pvec(a,b) ==> c = a x b
%
c(:,1)=a(:,2).*b(:,3)-a(:,3).*b(:,2);
c(:,2)=a(:,3).*b(:,1)-a(:,1).*b(:,3);
c(:,3)=a(:,1).*b(:,2)-a(:,2).*b(:,1);