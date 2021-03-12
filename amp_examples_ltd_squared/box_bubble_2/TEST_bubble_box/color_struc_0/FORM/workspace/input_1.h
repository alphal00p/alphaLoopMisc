CTable uvdiag(0:0);


L CONF =
 +conf(0,0,cmb(k1,+c4,k2,+c5,k3,+c2,k4,+c3,k5,-c1+(p1+p2)-c2-c3),c4,c5,1)*diag(0,0,k1,k2)*diag(0,1)
;

L F = 1*(1)*(1-p1.p1*xT^-1-2*p1.p1*xS^-1-2*p1.k1*xT^-1-4*p1.k1*xS^-1+2*p1.k3*xT^-1+2*p1.k3*xS^-1+2*p1.k4*xS^-1-2*k1.k1*xT^-1-2*k1.k1*xS^-1+2*k1.k3*xT^-1+2*k1.k3*xS^-1+2*k1.k4*xS^-1-k3.k3*xT^-1-k3.k3*xS^-1-2*k3.k4*xS^-1-k4.k4*xS^-1)
;