function [N,dN_dpsi,dN_deta,dN_dceda] = ffHex8 (psi,eta,ceda) 
  %-------------------------------------------------
  % funcion que calcula las funciones de forma de un elem hexahedrico trilineal de 8 nodos
  % si, eta, ceda son las coordenadas naturales, coord definidas dentro del elem master [-1:x:1]
  % para checkiar si est'a OK, hacer ffHex8(nodo) tiene q ser la prop delta
  % la verificacion que no falla p elem standar, sum(N)=1
  %-------------------------------------------------
  
  xnod = [-1 -1 -1;
	  1 -1 -1;
	  1 1 -1;
	  -1 1 -1;
	  -1 -1 1;
	  1 -1 1;
	  1 1 1 ;
	  -1 1 1];
  N = 1/8 * (1+xnod(:,1)*psi).*(1+xnod(:,2)*eta).*(1+xnod(:,3)*ceda);
  
  dN_dpsi  = 1/8 * (xnod(:,1)).*(1+xnod(:,2)*eta).*(1+xnod(:,3)*ceda);
  dN_deta  = 1/8 * (1+xnod(:,1)*psi).*xnod(:,2).*(1+xnod(:,3)*ceda);
  dN_dceda = 1/8 * (1+xnod(:,1)*psi).*(1+xnod(:,2)*eta).*xnod(:,3);
  
end
