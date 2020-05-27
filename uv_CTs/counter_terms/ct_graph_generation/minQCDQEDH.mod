* PDGS
*d,u,c,s,b,t -> 1 to 6
*d~,u~,c~,s~,b~,t~ -> -1 to -6
*g, photon -> 21, 22
*e+ e- > -11, 11
*mu+, mu-, ta+ ,ta- > -12, 12, -12, 13
*w+ w- z -> 24, -24, 23
*h -> 25 
* ghost -> 82

* propagators: bosons
[higgs, higgs, +; pdg=('25')]
[gluon,gluon,+; pdg=('21')]
[photon,photon,+;pdg=('22')]

* propagators: fermions, Psi PsiBar
* u-quark is massless uptype representative
* t-quark is massive  uptype representative

[ghost,ghostbar,-; pdg=('82','-82')]
[eminus, eplus, - , external; pdg= ('+11', '-11')]
[u, ubar, - ; pdg= ('+2','-2')]
[t, tbar, - ; pdg=('6','-6')]

* vertices qed 
[ubar,photon,u]
[tbar,photon,t]
[eplus,photon,eminus]

* vertices qcd + higgs
[ghostbar,gluon,ghost]
[ubar,gluon,u]
[tbar,gluon,t]
[gluon,gluon,gluon]
[gluon,gluon,gluon,gluon]
[tbar,higgs,t]
[higgs,higgs,higgs]
