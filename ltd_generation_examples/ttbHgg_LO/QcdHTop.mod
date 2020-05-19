*d,u,c,s,b,t -> 1 to 6
*d~,u~,c~,s~,b~,t~ -> -1 to -6
*g, photon -> 21, 22
*e+ e- > -11, 11
*mu+, mu-, ta+ ,ta- > -12, 12, -12, 13
*w+ w- z -> 24, -24, 23
*h -> 25
* propagators
[higgs,higgs,+;m='mh',pdg='25']
[gluon,gluon,+;m='0',pdg='21']
[u,ubar,-;m='0',pdg='2']
[t,tbar,-;m='mt',pdg='6']
[ghost,ghostbar,-;m='0',pdg='82']
[photon,photon,+;m='0',pdg='22']
[eminus,eplus,-,external;m='0',pdg='11']
*vertices
[ghostbar,gluon,ghost]
[ubar,gluon,u]
[tbar,gluon,t]
[gluon,gluon,gluon]
[gluon,gluon,gluon,gluon]
*[ubar,phot,u]
*[tbar,phot,t]
*[phot,phot,higgs]
[tbar,higgs,t]
[eplus,photon,eminus]
[tbar,photon,t]
*[higgs,gluon,gluon]
*[higgs,gluon,gluon,gluon]
*[higgs,gluon,gluon,gluon,gluon]
