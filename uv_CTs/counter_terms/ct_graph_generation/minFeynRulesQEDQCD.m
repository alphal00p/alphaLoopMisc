(* ::Package:: *)

(* Created with the Wolfram Language : www.wolfram.com *)
(* only numerator information, denominators are clear from particleType + momentumMap in diagram *)
SetAttributes[v,Orderless];
   
   prop[gluon[ p_,{n_,m_}]] := -(ii*g[mu[n], mu[m]]*delta[sunA[n],sunA[m]]) scalarProp[p,0]; 
   prop[photon[ p_,{n_,m_}]] := -(ii*g[mu[n], mu[m]])scalarProp[p,0]; 
   prop[u[ p_,{n_,m_}]] := ii*delta[sunF[n], sunF[m]]*gamma[s[n], mu[n, m], s[m]]*vector[p, mu[n, m]]scalarProp[p,0];
   prop[t[ p_,{n_,m_}]] := ii*delta[sunF[n], sunF[m]]*(gamma[s[n], mu[n, m], s[m]]*vector[p, mu[n, m]]+mt delta[s[n],s[m]])scalarProp[p,0];
   prop[ghost[p_,{n_,m_}]] := ii*delta[sunA[n], sunA[m]]scalarProp[p,0];
   prop[higgs[ p_,{n_,m_}]] := -ii;
   (* gluon-ghost *)
   v[ghost[p2_,i2_], ghostbar[p1_,i1_], gluon[p3_,i3_]] := gs*f[sunA[i3], sunA[i1], sunA[i2]]*(vector[p2, mu[i3]] + vector[p3, mu[i3]]);
   (* pure gluon *)
   v[gluon[p1_,i1_], gluon[p2_,i2_], gluon[p3_,i3_]] := (gs*f[sunA[i1], sunA[i2], sunA[i3]]*(-(g[mu[i1], mu[i3]]*vector[p1, mu[i2]]) + 
      g[mu[i1], mu[i2]]*vector[p1, mu[i3]] + g[mu[i2], mu[i3]]*
      vector[p2, mu[i1]] - g[mu[i1], mu[i2]]*vector[p2, mu[i3]] - 
      g[mu[i2], mu[i3]]*vector[p3, mu[i1]] + g[mu[i1], mu[i3]]*
      vector[p3, mu[i2]])); 
   v[gluon[p1_,i1_], gluon[p2_,i2_], gluon[p3_,i3_], gluon[p4_,i4_]] :> gs^2*ii*(f[sunA[i1], sunA[i2], sunA[i1, i2, i3, i4]]*f[sunA[i3], sunA[i4], sunA[i1, i2, i3, i4]] *(g[mu[i1], mu[i4]]*g[mu[i2], mu[i3]] - 
      g[mu[i1], mu[i3]]*g[mu[i2], mu[i4]]) + 
      f[sunA[i1], sunA[i3], sunA[i1, i2, i3, i4]]*f[sunA[i2], sunA[i4], sunA[i1, i2, i3, i4]]*
      (g[mu[i1], mu[i4]]*g[mu[i2], mu[i3]] - g[mu[i1], mu[i2]]*
      g[mu[i3], mu[i4]]) + f[sunA[i1], sunA[i4], sunA[i1, i2, i3, i4]]*
      f[sunA[i2], sunA[i3], sunA[i1, i2, i3, i4]]*
      (g[mu[i1], mu[i3]]*g[mu[i2], mu[i4]] - g[mu[i1], mu[i2]]*
       g[mu[i3], mu[i4]]));
   (* g q qbar *) 
   v[gluon[p3_,i3_], u[p2_,i2_], ubar[p1_,i1_]] := gs*ii*gamma[s[i1], mu[i3], s[i2]]*T[sunA[i3], sunF[i1], sunF[i2]];
   v[gluon[p3_,i3_], t[p2_,i2_], tbar[p1_,i1_]] := gs*ii*gamma[s[i1], mu[i3], s[i2]]*T[sunA[i3], sunF[i1], sunF[i2]];
   (* gamma q qbar *) 
   v[photon[p3_,i3_], u[p2_,i2_], ubar[p1_,i1_]] := charge[u] ee*ii*gamma[s[i1], mu[i3], s[i2]] delta[sunF[i1],sunF[i2]]; (* want ddB *)
   v[photon[p3_,i3_], t[p2_,i2_], tbar[p1_,i1_]] := charge[u] ee*ii*gamma[s[i1], mu[i3], s[i2]] delta[sunF[i1],sunF[i2]]; (* want ddB *)
   (*eplus eminus gamma*)      
   v[photon[p3_,i3_], eminus[p2_,i2_], eplus[p1_,i1_]] := (-1)*ee*ii*gamma[s[i1], mu[i3], s[i2]] ;
   (* t tbar higgs *)
   v[higgs[p2_,i2_],t[p1_,i1_],tbar[p3_,i3_]]:=- ii (yt/Sqrt[2]) delta[s[i1],s[i3]] delta[sunF[i1],sunF[i3]] 
   out[x__]:=1;
   in[x__]:=1;
	(*
   (* initial states *)
   in[eminus[p_,i_]]:=spinorU[{p,0},s[i]]
   in[eplus[p_,i_]]:=spinorVbar[{p,0},s[i]]
   (* final states *)
   out[eminus[p_,i_]]:=spinorUbar[{p,0},s[i]]
   out[eplus[p_,i_]]:=spinorV[{p,0},s[i]]

   (* initial states *)
   in[q[p_,i_]]:=spinorU[{p,0},s[i]]
   in[qbar[p_,i_]]:=spinorVbar[{p,0},s[i]]
   (* final states *)
   out[q[p_,i_]]:=spinorUbar[{p,0},s[i]]
   out[qpar[p_,i_]]:=spinorV[{p,0},s[i]]
*)
