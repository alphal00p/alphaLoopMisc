CTable uvdiag(0:0);


L CONF =
 +conf(0,0,cmb(k1,+c2,k2,+c3,k3,+c1),c2,c3,1)*diag(0,0,k1,k2)*diag(0,1)
;

L F = 1*(1)*(1-p1.p1*xT^-1-p1.p1*xS^-1-2*p1.k1*xT^-1-2*p1.k1*xS^-1+2*p1.k3*xT^-1-p2.p2*xS^-1+2*p2.k1*xS^-1-2*k1.k1*xT^-1-2*k1.k1*xS^-1+2*k1.k3*xT^-1-k3.k3*xT^-1)
;