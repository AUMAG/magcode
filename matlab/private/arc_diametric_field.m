function magB = arc_diametric_field(mag,g_xyz)

magB = [0;0;0];

thetaS = 0;
theta  = 0;
rho    = 0;


for mm = 1:2
  for nn = 1:2
    for qq = 1:2
      magB = (-1)*(mm+nn+qq)*[ ...
        Brho1(mm,nn,qq)+Brho2(mm,nn,qq);...
        Bthe1(mm,nn,qq)+Bthe2(mm,nn,qq);...
        Bzzz1(mm,nn,qq)+Bzzz2(mm,nn,qq)...
      ];
    end
  end
end

magB = magB*MT*(2e7);


function c = Brho1(mm,nn,qq)

  rhod_m = 0;
  varrho_m = 0;
  barrho_m = 0;
  Z_n = 0;
  Ghi_mnq = 0;
  R_mn = 0;
  kkappa_m = 0;
  phi_q = 0;
  kk_nm = 0;

  aa = -(1+2*barrho_m/(varrho_m*kkappa_m))*EllipticF(phi_q,kk_nm);
  bb = 2*EllipticD(phi_q,kk_nm);
  cc = barrho*(kkappa_m-2)/(varrho_m*kkappa_m)*EllipticPi(phi_q,kkappa_m,kk_nm);

  a = cos(thetaS-theta)*rhod_m*Z_n/(rho*R_mn)*(aa+bb+cc);
  b = sin(thetaS-theta)/(2*rho^2)*(varrho_m*barrho_m*atanh(Ghi_mnq/Z_n)-Z_n*Ghi_mnq);
  c = a+b;

end

end
