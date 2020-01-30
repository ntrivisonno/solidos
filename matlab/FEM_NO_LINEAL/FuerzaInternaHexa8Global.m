function F = FuerzaInternaHexa8Global (XYZ,CONEC,U,funcTension)
	 
  nnod = size(XYZ,1);
  nelem = size(CONEC,1);
  
  F = zeros(3*nnod,1);
  % grados de libertad del elem
  glo = zeros(1,24);
  
  for e = 1:nelem
    % nodos del elem
    conec = CONEC(e,:);
    for j = 1:3
      glo(j:3:end) = 3*(conec-1)+j;
    end
    % coord globales de los nodos del elem
    x = XYZ(conec,1);
    y = XYZ(conec,2);
    z = XYZ(conec,3);
    
    % desplazamiento de los nodos del elemento
    u = U(glo);
    
    % calc de fuerza interna
    f = FuerzaInternaHexa8 (x,y,z,u,funcTension); 
    F(glo) = F(glo) + f; 
  end

end
