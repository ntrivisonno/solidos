function K = TangenteHexa8Global (XYZ,CONEC,U,funcTension,funcMod,calcid)
  % funcMod -> C	 

  nnod = size(XYZ,1);
  nelem = size(CONEC,1);
  
  K = zeros(3*nnod);
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
    
    % desplazamientos de los nodos del elemento
    u = U(glo);
    
    % fuerzas internas en el elemento
    if calcid ==1 % calc tangente numerica
      k = TangenteNumHexa8(x,y,z,u,funcTension);
    elseif calcid == 2 % tangente analitica
      k = TangenteHexa8(x,y,z,u,funcTension,funcMod);
    elseif calcid ==3 % rigidez lineal
      k = rigidezLinealHexa8 (x,y,z,funcMod);
    end
    % ensamblaje
    K(glo,glo) = K(glo,glo) + k; 
  end

end
