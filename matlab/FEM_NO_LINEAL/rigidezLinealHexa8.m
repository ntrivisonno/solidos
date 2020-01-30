function k = rigidezLinealHexa8 (x,y,z,C)
  % C -> tensor de mod elastico
	 
  npg = 8;
  nnodxelem = 8;
  ndofxelem = nnodxelem*3;
  k = zeros(ndofxelem);
  for i = 1:npg
    % peso de coord naturales del punto de integracion
    [w, psi, eta, ceda] = WeightAndCoordHexa(npg,i);
    [J,dN_dpsi,dN_deta,dN_dceda] = jacobiano (x,y,z,psi,eta,ceda);
    j = inv(J);
    dNdx = zeros ( nnodxelem,1);
    dNdy = zeros ( nnodxelem,1);
    dNdz = zeros ( nnodxelem,1);
    for k1 = 1:nnodxelem
      dNknat = [dN_dpsi(k1);dN_deta(k1);dN_dceda(k1)];
      dNkglob = j*dNknat;
      dNdx(k1) = dNkglob(1);
      dNdy(k1) = dNkglob(2);
      dNdz(k1) = dNkglob(3);
    end
    
    Blin = zeros(6,24);
    Blin(1,1:3:end) = dNdx;
    Blin(2,2:3:end) = dNdy;
    Blin(3,3:3:end) = dNdz;
    Blin(4,1:3:end) = dNdy;
    Blin(4,2:3:end) = dNdx;
    Blin(5,2:3:end) = dNdz;
    Blin(5,3:3:end) = dNdy;
    Blin(6,1:3:end) = dNdz;
    Blin(6,3:3:end) = dNdx;
    
    f = Blin'*C*Blin;
    k = k + f*w*det(J);
  end
  
end
