function mtl2vtu(FILE,NOD,ELE,SOL)
  % mtl2vtu(F2WRITE,NOD,ELE,SOL)
  % Writes file F2WRITE in vtu format.
  % NOD can be a structure or a matrix with nodal coordinates
  %    a) if ismatrix(NOD), NOD(i,j): j-coordinate of node i (j=1,2(,3), i=1,2,...,nnod)
  %    b) is isstruct(NOD), NOD must have the field Coordinates.
  %                         NOD.Coordinates(i,j): j-coordinate of node i
  % ELE is a structure defining the element set, which must have
  %        have the fields Connectivity and Type:
  %                         ELE.Connectivity(i,j): j-th node of the element i
  %                         ELE.Type: type of finite element
  % SOL(i) is a structure (i=1,...,nbr of quantities to display) with 4 fields:
  %                   SOL(i).Name
  %                   SOL(i).Type (PointData or CellData)
  %                   SOL(i).NumberOfComponents
  %                   SOL(i).Value
  % If SOL is empty or does not exist, only the mesh is displayed.
  %
	 
  FID = fopen(FILE,'w');
  fprintf(FID, '<?xml version="1.0"?>\n');
  fprintf(FID, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
  fprintf(FID, '    <UnstructuredGrid>\n');
  
  if isstruct(NOD)
    NOD = NOD.Coordinates;
  end
  
  [nnod,dim] = size(NOD);
  if dim==2
    NOD(1,3) = 0;
    
  end
  NOD = NOD';
  
  ETYPE = ELE.Type;
  CONEC = ELE.Connectivity;
  [nele,nnodperele] = size(CONEC);
  fprintf(FID, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n',nnod,nele);
  fprintf(FID, '<Points>\n');
  fprintf(FID, '<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
  fprintf(FID, '%f ',NOD(:)');
  fprintf(FID, '\n</DataArray>\n');
  fprintf(FID, '</Points>\n');
  
  switch(ETYPE)
    case {'BAR2'}
      etype = 3; %VTK_LINE
    case {'TRI3'}
      etype = 5; %VTK_TRIANGLE
    case {'TRI6'}
      etype = 22; %VTK_QUADRATIC_TRIANGLE
    case {'QUAD4','MITC4'}
      etype = 9; % VTK_QUAD
    case {'QUAD8'}
      etype = 23; %VTK_QUADRATIC_QUAD
    case {'QUAD9'}
      etype = 28; %VTK_BIQUADRATIC_QUAD
    case {'QUAD16'}
      % etype = 36; %VTK_BICUBIC_QUAD
      % No funciona. Paraview 3.12.0 no acepta cuadrángulos bi-cúbicos.
      % Para visualización, cada cuadrángulo bi-cúbico es dividido en
      % 9 cuadrángulos bi-lineales.
      etype = 9; % VTK_QUAD
      C = CONEC;
      CONEC = [C(:,[1 5 13 12])
               C(:,[5 6 14 13])
               C(:,[6 2 7 14])
               C(:,[12 13 16 11])
               C(:,[13 14 15 16])
               C(:,[14 7 8 15])
               C(:,[11 16 10 4])
               C(:,[16 15 9 10])
               C(:,[15 8 3 9])];
      nele = 9*nele;
      nnodperele = 4;
    case {'TETRA4'}
      etype = 10; %VTK_TETRA
    case {'TETRA10'}
      % vtkQuadraticTetra is a concrete implementation of vtkNonLinearCell to represent a three-dimensional, 10-node, isoparametric parabolic tetrahedron.
      % The interpolation is the standard finite element, quadratic isoparametric shape function.
      % The cell includes a mid-edge node on each of the size edges of the tetrahedron. The ordering of the ten points defining the cell is
      % point ids (0-3,4-9) where ids 0-3 are the four tetra vertices; and point ids 4-9 are the midedge nodes between (0,1), (1,2), (2,0), (0,3), (1,3), and (2,3).
      etype = 24;
      CONEC = CONEC(:,[1 2 3 4 5 8 6 7 10 9]);
    case {'HEXA8','HEXA8RI','HEXA8MITC4','VOLHEXA8'}
      etype = 12; %VTK_HEXAHEDRON
    case {'HEXA20'}
      etype = 25; %VTK_QUADRATIC_HEXAHEDRON
    case {'HEXA27'}
      etype = 29; %VTK_TRIQUADRATIC_HEXAHEDRON
      CONEC = CONEC(:,[1:20 24 26 23 25 21 22 27]);
    %             nnodperele = 20;
    %             ELE = ELE(:,1:20);
    otherwise
      error('Unknown element')
  end
  
  CONEC = CONEC';
  offs = nnodperele:nnodperele:nnodperele*nele;
  typs = etype*ones(1,nele);

  fprintf(FID, '<Cells>\n');
  fprintf(FID, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
  fprintf(FID,'%d ',CONEC(:)'-1);
  fprintf(FID,'\n</DataArray>\n');
  fprintf(FID,'<DataArray type="Int32" Name="offsets" format="ascii">\n');
  N = num2str(offs,' %d');
  fprintf(FID,'%d ',offs);
  fprintf(FID, '\n</DataArray>\n');
  fprintf(FID, '<DataArray type="UInt8" Name="types" format="ascii">\n');
  fprintf(FID, '%d ',typs);
  fprintf(FID, '\n</DataArray>\n');
  fprintf(FID, '</Cells>\n');


  if nargin==4 && isstruct(SOL)
    nsol = length(SOL);
    pdata = [];
    cdata = [];
    for i=1:nsol
      if strcmp(SOL(i).Type,'PointData')
        pdata = [pdata i];
        else
          cdata = [cdata i];
      end
    end
    if ~isempty(pdata)
      fprintf(FID, '<PointData>\n');
      for i = pdata
        sol = SOL(i);
        fprintf(FID,'<DataArray type="Float32" Name="%s " NumberOfComponents="%d " format="ascii">\n',sol.Name,sol.NumberOfComponents);
        v = sol.Value;
        if issparse(v)
          v = full(v);
        end
        [n,d] = size(v);
        if d > 1
          if sol.NumberOfComponents ~= d
            error('Wrong number of components')
          end
          V = reshape(v',1,[]);
        else
          V = v(:)';
        end
        fprintf(FID,'%f ',V);
        fprintf(FID,'\n</DataArray>\n');
      end
      fprintf(FID, '</PointData>\n');
    end
    
    
    if ~isempty(cdata)
      fprintf(FID, '<CellData>\n');
      for i = cdata
        sol = SOL(i);
        fprintf(FID,'<DataArray type="Float32" Name="%s " NumberOfComponents="%d " format="ascii">\n',sol.Name,sol.NumberOfComponents);
        v = sol.Value;
        if issparse(v)
          v = full(v);
        end
        [n,d] = size(v);
        if d > 1
          if sol.NumberOfComponents ~= d
            error('Wrong number of components')
          end
          V = reshape(v',1,[]);
        else
          V = v(:)';
        end
        fprintf(FID,'%f ',V);
        fprintf(FID,'\n</DataArray>\n');
      end
      fprintf(FID, '            </CellData>\n');
    end
  end
  
  fprintf(FID, '        </Piece>\n');
  fprintf(FID, '    </UnstructuredGrid>\n');
  fprintf(FID, '</VTKFile>\n');
  
  fclose(FID);
%toc
end
