(* The finite subtracted one-loop integrand. The numerical value of p1 must be proportional to (1,0,0,1), and the numerical value of p2 must be proportional to (1,0,0,-1).
p1 and p2 are the incoming electron / positron momenta. p3 and p4 are the outgoing off-shell photon momenta, and the corresponding photon polarization vectors are written as p5 and p6. We must have p1+p2=p3+p4. k1 is the loop momentum.*)
{(-I)*d[mu1, mu2]*d[mu13, mu14]*gamma[tree, p2, mu14, k1 - p2, p6, 
   k1 + p1 - p3, p5, k1 + p1, mu13, p1]*propagator[k1]*propagator[k1 + p1]*
  propagator[k1 - p2]*propagator[k1 + p1 - p3], 
 (I*d[mu1, mu2]*d[mu13, mu14]*gamma[tree, p2, mu14, k1 - p2, p1, p2, p6, 
    p1 - p3, p5, p1, p2, k1 + p1, mu13, p1]*propagator[k1]*
   propagator[k1 + p1]*propagator[k1 - p2]*propagator[p1 - p3])/
  lorentzDot[p1 + p2]^2, (-I)*d[mu1, mu2]*d[mu13, mu14]*
  gamma[tree, p2, p6, p1 - p3, mu14, k1 + p1 - p3, p5, k1 + p1, mu13, p1]*
  propagator[k1]*propagator[k1 + p1]*propagator[p1 - p3]*
  propagator[k1 + p1 - p3], (-I)*d[mu1, mu2]*d[mu11, mu12]*
  gamma[tree, p2, p6, p1 - p3, mu12, k1 + p1 - p3, mu11, p1 - p3, p5, p1]*
  propagator[k1]*propagator[p1 - p3]^2*propagator[k1 + p1 - p3], 
 (-I)*d[mu1, mu2]*d[mu13, mu14]*gamma[tree, p2, mu13, k1 - p2, p6, 
   k1 + p1 - p3, mu14, p1 - p3, p5, p1]*propagator[-k1]*propagator[k1 - p2]*
  propagator[p1 - p3]*propagator[k1 + p1 - p3], 
 (-I)*d[mu1, mu2]*d[mu12, mu11]*gamma[tree, p2, p6, p1 - p3, mu12, k1, 
   p1 - p3, k1, mu11, p1 - p3, p5, p1]*propagator[p1 - p3]^2*
  propagator[k1, M]^3, I*d[mu1, mu2]*d[mu12, mu11]*
  gamma[tree, p2, p6, p1 - p3, mu12, k1, mu11, p1 - p3, p5, p1]*
  propagator[p1 - p3]^2*propagator[k1, M]^2, I*d[mu1, mu2]*d[mu13, mu14]*
  gamma[tree, p2, mu13, -k1, p6, -k1, mu14, p1 - p3, p5, p1]*
  propagator[p1 - p3]*propagator[-k1, M]^3, 
 (I*d[mu1, mu2]*d[mu14, mu13]*
   (-gamma[tree, p2, mu14, k1, p1, p2, p6, p1 - p3, p5, p1, p2, k1, mu13, 
      p1] + gamma[tree, p2, p6, p1 - p3, mu14, k1, p5, k1, mu13, p1]*
     lorentzDot[p1 + p2]^2)*propagator[p1 - p3]*propagator[k1, M]^3)/
  lorentzDot[p1 + p2]^2}
