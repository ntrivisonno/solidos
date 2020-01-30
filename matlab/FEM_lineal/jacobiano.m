function [J,dN_dpsi,dN_deta,dN_dceda] = jacobiano (x,y,z,psi,eta,ceda)
  %-------------------------------------------------
  % x y z son las coordenadas de nodo global 
  % son las coord de nodo natural
  %-------------------------------------------------
	 
  [~,dN_dpsi,dN_deta,dN_dceda] = ffHex8 (psi,eta,ceda);
  J = zeros(3);
  J = [dot(dN_dpsi,x)   dot(dN_dpsi,y)   dot(dN_dpsi,z);
	    dot(dN_deta,x)   dot(dN_deta,y)   dot(dN_deta,z);
	    dot(dN_dceda,x)  dot(dN_dceda,y)  dot(dN_dceda,z)];
  
end
