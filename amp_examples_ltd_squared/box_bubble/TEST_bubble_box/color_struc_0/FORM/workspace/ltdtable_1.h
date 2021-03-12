
Auto S invd, E, shift, ltd;
S r, s;
CF a, num, ncmd, conf1, replace, energies, ellipsoids, ltdcbtolmb, ltdenergy, constants;
NF allenergies;
Set invdset: invd0,...,invd10000;
CTable ltdtopo(0:0);

Fill ltdtopo(0) = (-1)^2*constants()*
	ellipsoids(2*E0,+E0-energies(0)+E1+energies(p1),+E0-energies(0)-E1+energies(p1),+E0-energies(0)+E2+energies(-k3+p1),+E0-energies(0)-E2+energies(-k3+p1),+E0-energies(0)+E3+energies(-p2),+E0-energies(0)-E3+energies(-p2),2*E4,+E1-energies(p1)+E0+energies(0),+E1-energies(p1)-E0+energies(0),2*E1,+E1-energies(p1)+E2+energies(-k3+p1),+E1-energies(p1)-E2+energies(-k3+p1),+E1-energies(p1)+E3+energies(-p2),+E1-energies(p1)-E3+energies(-p2),2*E4,+E2-energies(-k3+p1)+E0+energies(0),+E2-energies(-k3+p1)-E0+energies(0),+E2-energies(-k3+p1)+E1+energies(p1),+E2-energies(-k3+p1)-E1+energies(p1),2*E2,+E2-energies(-k3+p1)+E3+energies(-p2),+E2-energies(-k3+p1)-E3+energies(-p2),2*E4,+E3-energies(-p2)+E0+energies(0),+E3-energies(-p2)-E0+energies(0),+E3-energies(-p2)+E1+energies(p1),+E3-energies(-p2)-E1+energies(p1),+E3-energies(-p2)+E2+energies(-k3+p1),+E3-energies(-p2)-E2+energies(-k3+p1),2*E3,2*E4)*
	allenergies(k1,k1+p1,k1-k3+p1,k1-p2,k2)*(
		(1*ltdcbtolmb(c2,ltd0-energies(0),c3,ltd4-energies(0))*ltdenergy(ltd0,+E0,ltd4,+E4)*der(ltd4,1)*(prop(1,1,ltd0,0,E0)*prop(2,1,ltd0,-energies(0),E1+energies(p1))*prop(3,1,ltd0,-energies(0),-E1+energies(p1))*prop(4,1,ltd0,-energies(0),E2+energies(-k3+p1))*prop(5,1,ltd0,-energies(0),-E2+energies(-k3+p1))*prop(6,1,ltd0,-energies(0),E3+energies(-p2))*prop(7,1,ltd0,-energies(0),-E3+energies(-p2))*prop(8,1,ltd4,0,E4)^2)
			+1*ltdcbtolmb(c2,ltd1-energies(p1),c3,ltd4-energies(0))*ltdenergy(ltd1,+E1,ltd4,+E4)*der(ltd4,1)*(prop(9,1,ltd1,-energies(p1),E0+energies(0))*prop(10,1,ltd1,-energies(p1),-E0+energies(0))*prop(11,1,ltd1,0,E1)*prop(12,1,ltd1,-energies(p1),E2+energies(-k3+p1))*prop(13,1,ltd1,-energies(p1),-E2+energies(-k3+p1))*prop(14,1,ltd1,-energies(p1),E3+energies(-p2))*prop(15,1,ltd1,-energies(p1),-E3+energies(-p2))*prop(16,1,ltd4,0,E4)^2)
			+1*ltdcbtolmb(c2,ltd2-energies(-k3+p1),c3,ltd4-energies(0))*ltdenergy(ltd2,+E2,ltd4,+E4)*der(ltd4,1)*(prop(17,1,ltd2,-energies(-k3+p1),E0+energies(0))*prop(18,1,ltd2,-energies(-k3+p1),-E0+energies(0))*prop(19,1,ltd2,-energies(-k3+p1),E1+energies(p1))*prop(20,1,ltd2,-energies(-k3+p1),-E1+energies(p1))*prop(21,1,ltd2,0,E2)*prop(22,1,ltd2,-energies(-k3+p1),E3+energies(-p2))*prop(23,1,ltd2,-energies(-k3+p1),-E3+energies(-p2))*prop(24,1,ltd4,0,E4)^2)
			+1*ltdcbtolmb(c2,ltd3-energies(-p2),c3,ltd4-energies(0))*ltdenergy(ltd3,+E3,ltd4,+E4)*der(ltd4,1)*(prop(25,1,ltd3,-energies(-p2),E0+energies(0))*prop(26,1,ltd3,-energies(-p2),-E0+energies(0))*prop(27,1,ltd3,-energies(-p2),E1+energies(p1))*prop(28,1,ltd3,-energies(-p2),-E1+energies(p1))*prop(29,1,ltd3,-energies(-p2),E2+energies(-k3+p1))*prop(30,1,ltd3,-energies(-p2),-E2+energies(-k3+p1))*prop(31,1,ltd3,0,E3)*prop(32,1,ltd4,0,E4)^2))
		*(1*(1))
);

