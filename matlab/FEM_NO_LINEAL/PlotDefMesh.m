function PlotDefMesh(coord,conec,u,sc,etype,lin)
  % PlotMesh(coord,conec,u,sc,etype,lin)
  % Function that plots a finite element mesh.
  % coord(i,j): j-coordinate of node i 
  % conec(i,j): j-th node of element i
  % u(i,j): j-displacement of node i 
  % sc: deformed scale
  % etype = 'BAR2', 'TRI3', ...
  % lin: line type (by default, lin = 'r-' (solid red line))
	 
  if nargin == 5
    lin = 'r-';
  end
  
  [nnod,dim] = size(coord);
  
  switch(etype)
    case {'BAR2', 'Bar2'}
      v = [1 2];
    case {'TRI3', 'Tri3'}
      v = [1 2 3 1];        
    case {'TRI6', 'Tri6'}
      v = [1 4 2 5 3 6 1];        
    case {'QUAD4', 'Quad4'}
      v = [1 2 3 4 1];
    case {'QUAD8', 'Quad8','QUAD9', 'Quad9'}
      v = [1 5 2 6 3 7 4 8 1];
    case {'QUAD16', 'Quad16'}
      v = [1 5 6 2 7 8 3 9 10 4 11 12 1];
    case {'TETRA4', 'Tetra4'}
      v = [1 2 4 1  % 2 3 4 2
           1 4 3 1];
    case {'HEXA8', 'Hexa8'}
      v = [1 2 3 4 1
           5 6 7 8 5
           1 2 6 5 1  % 2 3 7 6 2
           3 4 8 7 3  % 4 1 5 8 4
          ];    
    case {'HEXA20', 'Hexa20', 'HEXA27', 'Hexa27'}
      v = [1 9 2 10 3 11 4 12 1
           5 13 6 14 7 15 8 16 5
           1 9 2 18 6 13 5 17 1   % 2 10 3 19 7 14 6 18 2
           3 11 4 20 8 15 7 19 3  % 4 12 1 17 5 16 8 20 4
          ];    
    otherwise
      error('Not coded yet')
  end
  [nele,nnodperele] = size(conec);
  [nf,nv] = size(v);
  X = zeros(nv,nele);
  Y = X;
  Z = X;
  if size(u,2)==dim
    U = sc*u;
  else
    U = sc*reshape(u,dim,nnod)';
  end
  
  for j = 1:nf
    hold on
    for i = 1:nv
      V = conec(:,v(j,i));
      X(i,:) = (coord(V,1)+U(V,1))';
      Y(i,:) = (coord(V,2)+U(V,2))';
      if dim == 3
        Z(i,:) = (coord(V,3)+U(V,3))';
      end
    end
    if dim==2
      plot(X,Y,lin);
    else
      plot3(X,Y,Z,lin);
    end
    hold off
  end

end
