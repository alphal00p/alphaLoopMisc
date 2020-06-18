(* ::Package:: *)

removeNestedFermionChain[expr_] := expr /. id -> id[4] /. {
    fermionChain[a__] :> (Expand[#, gamma[__]|slash[__]|id[__]| spinor[__]|projectorP1P2]& /@ fermionChain[a])
    } //.{
    fermionChain[rest1___, a_/;FreeQ[a, gamma[__] | slash[__] |id[__]| spinor[__]|projectorP1P2 | fermionChain[__]], rest2___] :>
        a * fermionChain[rest1, rest2],
    fermionChain[rest1___, a_+b_ /; Not[FreeQ[a, fermionChain[__]]], rest2___] :>
        fermionChain[rest1, a, rest2] +
        fermionChain[rest1, b, rest2],
    fermionChain[rest1___, factor1_*fermionChain[b__], rest3___]:>
        factor1 * fermionChain[rest1, b, rest3],
    fermionChain[rest1___, fermionChain[b__], rest3___]:>
        fermionChain[rest1, b, rest3]
};


keepOnlyLoopMomenta[expr_] := expr /. {p1->0, p2->0, q1->0, q2->0} //. {lower[0] -> 0, minusPart[0]->0, plusPart[0] ->0};(*Module[{expr1},
    expr1 = expr /. {lower[k1]\[Rule]lowerK1, lower[k2]\[Rule]lowerK2};
    D[expr1, k1] * k1 + 
    D[expr1, k2] * k2 +
    D[expr1, lowerK1] * lowerK1 +
    D[expr1, lowerK2] * lowerK2 +
    D[expr1, upperK1] * upperK1 +
    D[expr1, upperK2] * upperK2
] /. {lowerK1\[Rule]lower[k1], lowerK2\[Rule]lower[k2]};*)


keepOnlyExternalMomenta[expr_]:=
    expr /. {k1->0, k2->0};


noLoopQ[mom_] := FreeQ[mom, k1 | k2];


insertMassRegulatorAndUVexpand[expr_] :=
    expr /. slash[a_] :> slash[Expand[a]] //. slash[b_+c_] :> slash[b] + slash[c] /.
        propagator[mom_] :> propagator[mom, 0] /. {propagator[mom_, mass_] :>
        If[noLoopQ[mom],
            propagator[mom, mass],
            propagator[keepOnlyLoopMomenta[mom], M] -
                propagator[keepOnlyLoopMomenta[mom], M]^2 * 2 * lorentzDot[keepOnlyLoopMomenta[mom], keepOnlyExternalMomenta[mom]]
        ]
    } /. {lorentzDot[_,0]->0, (lorentzDot|component)[0,_]->0};


breakUpFermionChain[expr_] :=
    expr /.
     {fermionChain[a__] :> Module[ {list},
                               list = List[a];
                               list = 
                                Table[list[[i]] /. {obj: (slash|gamma|id|spinor)[b__] | projectorP1P2 :> 
                                    diracGamma[i, obj]}, {i, Length[list]}];
                               Times @@ list
                           ]
      };


rebuildFermionChain[expr_] :=
    expr /. {diracGamma[i_, b1_] :> 
         fermionChain[diracGamma[i, b1]]} //.
      {
       fermionChain[rest1___, diracGamma[i_, b1_]] * 
          fermionChain[diracGamma[j_, b2_], rest2___] /; j == i + 1 :> 
        fermionChain[rest1, diracGamma[i, b1], diracGamma[j, b2], 
         rest2]
       } /.
     {diracGamma[i_Integer, a_] :> a};


uvScalingDiagram[expr_] :=
    expr /. lorentzDot[a_] :> lorentzDot[a,a] /. {propagator[mom_, mass_] :>
       If[noLoopQ[mom],
           propagator[mom, mass],
           propagator[mom, mass] * delta^2
       ],
      diracGamma[i_Integer, slash[a_]] :>
       If[ noLoopQ[a],
           diracGamma[i, slash[a]],
           1 / delta * diracGamma[i, slash[keepOnlyLoopMomenta[a]]]
       ],
      (head:lorentzDot|component)[a_, b_] :>
       If[ noLoopQ[a] && noLoopQ[b],
           head[a, b],
           If[noLoopQ[a],
               1 / delta * head[a, keepOnlyLoopMomenta[b]],
               If[noLoopQ[b],
                   1 / delta * head[keepOnlyLoopMomenta[a], b],
                   1 / delta^2 * head[keepOnlyLoopMomenta[a], keepOnlyLoopMomenta[b]]
               ]
           ]
       ]
      };
