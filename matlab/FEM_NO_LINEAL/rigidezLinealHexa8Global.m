function K = rigidezLinealHexa8Global (XYZ,CONEC,C)
  % C -> tensor de mod elastico
	 
  nnod = size(XYZ,1);
  nelem = size(CONEC,1);
  
  K = zeros(3*nnod);
  
  for e = 1:nelem
    % nodos del elem
    conec = CONEC(e,:);
    % coord globales de los nodos del elem
    x = XYZ(conec,1);
    y = XYZ(conec,2);
    z = XYZ(conec,3);
    k = rigidezLinealHexa8 (x,y,z,C); % matriz de rigidez elem
    glo = zeros(1,24);
    
    for j = 1:3
      glo(j:3:end) = 3*(conec-1)+j;
    end
    K(glo,glo) = K(glo,glo) + k;
  % volumen del elem
    
  end

end
