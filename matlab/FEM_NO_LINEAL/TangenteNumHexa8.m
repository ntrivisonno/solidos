function K = TangenteNumHexa8 (x,y,z,U,funcTension)
  % funcTension -> tension en funcion de las relaciones constitutivas
	 
  F = FuerzaInternaHexa8(x,y,z,U,funcTension);
  K=zeros(24);
  dU = 1e-6;
  for k = 1:24
    U(k) =  U(k) + dU;
    F1 = FuerzaInternaHexa8(x,y,z,U,funcTension);
    K(:,k) = (F1-F)/dU;
    U(k) =  U(k) - dU;
  end
  
end
