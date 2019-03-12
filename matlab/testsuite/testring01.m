
clear all

r11 = 0.04;
r12 = 0.05;
h1  = 0.04;

r21 = 0.025;
r22 = 0.035;
h2  = 0.045;

cylmag11 = magnetdefine('type','cylinder','radius',r11,'height',h1,'grade','N42','dir','-z');
cylmag12 = magnetdefine('type','cylinder','radius',r12,'height',h1,'grade','N42','dir','+z');

cylmag21 = magnetdefine('type','cylinder','radius',r21,'height',h2,'grade','N42','dir','-z');
cylmag22 = magnetdefine('type','cylinder','radius',r22,'height',h2,'grade','N42','dir','+z');

ringmag1 = magnetdefine('type','cylinder','radius',[r11 r12],'height',h1,'grade','N42','dir','+z');
ringmag2 = magnetdefine('type','cylinder','radius',[r21 r22],'height',h2,'grade','N42','dir','+z');


%% Equal radii & height

N = 20;
zdispl = 0.1*(2*rand(N,1)-1);
zdispl = [0; h1/2+h2/2; zdispl];
 
for ii = 1:N
  
  F11 = magnetforces(cylmag11,cylmag21,[0; 0; zdispl(ii)]);
  F12 = magnetforces(cylmag11,cylmag22,[0; 0; zdispl(ii)]);
  F21 = magnetforces(cylmag12,cylmag21,[0; 0; zdispl(ii)]);
  F22 = magnetforces(cylmag12,cylmag22,[0; 0; zdispl(ii)]);
  
  Fcalc = F11 + F12 + F21 + F22;
  Fring = magnetforces(ringmag1,ringmag2,[0; 0; zdispl(ii)]);
  
  assert(all(Fcalc==Fring),['ring magnet calculations are consistent with cylindrical (z = ',num2str(zdispl(ii)),')'])

end