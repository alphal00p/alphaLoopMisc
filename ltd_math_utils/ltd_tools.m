(* ::Package:: *)

(* ::Input::Initialization:: *)
Print["The documented functions in this package are: \n ?contractLorentz \n ?contractColor \n Make sure you have the \"trace_form.log\" file for the gamma-algebra \n "];

contractLorentz::usage="contractLorentz[expr, optionalStatements]: 
	contractLorentz: Contracts Lorentz-indices and Spinor indices and performs Gamma-traces. It works with the following objects:
		Pure Lorentz Algebra:
		vector[p, mu]: defines the Lorentz vector \!\(\*SuperscriptBox[\(p\), \(\[Mu]\)]\)
		SP[p,k]: is the scalar-product \!\(\*SubscriptBox[\(p\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\(k\), \(\[Mu]\)]\)
		g[mu,nu]: is the metric tensor \!\(\*SubscriptBox[\(g\), \(\[Mu], \[Nu]\)]\)
		------------------------------------------------------------------
		Dirac Algebra:
		delta[s1,s2]: Identity in spinor space
		gamma[s1,mu,s2]: is \!\(\*SubscriptBox[\((\*SuperscriptBox[\(\[Gamma]\), \(\[Mu]\)])\), \(s1, s2\)]\) where s1,s2 are spinor indices
		gamma[s1,a,b,...,s2]: Gamma-chain where s1,s2 are the non-contracted spinor indices
		gamma[s1,lVec[p],s2]: is pSlash
		gamma[s1,mu,lVec[p],...,s2]: Gamma-chain where s1,s2 are the non-contracted spinor indices, slashed momenta are represented by lVec[momentum]
		-------------------------------------------------------------------
		Spinors:
		spinorUbar[{momentum,mass,polarization},spinorIndex]
		spinorVbar[{momentum,mass,polarization},spinorIndex]
		spinorU[{momentum,mass,polarization},spinorIndex]
		spinorV[{momentum,mass,polarization},spinorIndex]
	INPUT:
		expr: The expression to be contracted
		Optional Statements: (see: Options[contractLorentz] )
		- additionalRules\[Rule]{Userdefined rules e.g. SP[p1,p1]->0}
		- spinChainSimplify -> Boolean (default: False): If set to true, gamma-chains are brought into canonical (minimal) ordering, Dirac-equation is applied if possible
		- polarizationSum   -> Boolean ( default: False): If set to true, the polarization sum on external spinors is performed.		

	OUTPUT: Expression with contracted Lorentz indices.
	NEW CONSTANTS in Output: d for the dimension.
";
contractColor::usage=" contractColor[expr, quarkIndex , gluonIndex ]
	contractColor does the color-algebra. It works with the following objects:
		delta[gluonIndex[1],gluonIndex[2]],delta[quarkIndex[1],quarkIndex[2]] : Kronecker-delta k,n either \"quark-\" or \"gluon\"-indices
		f[gluonIndex[1],gluonIndex[2],gluonIndex[3]]: structure constant of su(Nc) (a,b,c): \"gluon-indices\"
		T[a,l,m]: (l,m)-th entry of the a-th generator of fundamental representation with (a): \"gluon-index \" and components (l,m): \"quark-indices\"
	INPUT: 
		expr: The expression on which color algebra should be performed
		quarkIndex: Symbol which is used to denote the entries of the generators of the fund. representation
		gluonIndex: Symbol used  to denote the a-th generator of the fundamental/adjoint-representation
					E.g. quarkIndex=i, gluonIndex=a expects generators like T[a[1],i[2],i[3]] and f[a[1],a[3],a[9]]
	OUTPUT: Expression with performed color-algebra.
	NEW CONSTANTS: Nc from su(Nc), TF from normalization of su(Nc) (e.g. TF=1/2)
";


Protect[additionalRules,traceFileLocation,spinChainSimplify,polarizationSum];
Options[contractLorentz] = {additionalRules -> {1->1}, traceFileLocation->$fileLocation,spinChainSimplify->False,polarizationSum->False }


PrintTemporary["load gamma traces"]

$fileLocation=FileNameJoin[{DirectoryName @ $InputFileName, ""}];

FromFormToMathe[filename_]:=Module[{imp},
		imp=StringJoin[ReadList[filename,String]];
		imp=StringReplace[imp,{"\\"->"","\n"->"","\r" -> "", "("->"[", ")"->"]"," "->"" ,"d_"->"g","mu"->"traceInd"}];
		ToExpression[StringSplit[imp,";"]]
	];





spinorSum[expr_]:=Block[{mom,mass,delta,pol,dummyLong,dummyLong2,lVec,sInd1,sInd2,res},
    gamma[spinorUbar[{mom_,mass_,pol_}],dummyLong___,spinorU[{mom_,mass_,pol_}]]=gammaTrace[lVec[mom],dummyLong]+mass gammaTrace[dummyLong];
    gamma[spinorVbar[{mom_,mass_,pol_}],dummyLong___,spinorV[{mom_,mass_,pol_}]]=gammaTrace[lVec[mom],dummyLong]-mass gammaTrace[dummyLong];
    
    gamma[dummyLong__,spinorU[{mom_,mass_,pol_}]]gamma[spinorUbar[{mom_,mass_,pol_}],dummyLong2__]^=gamma[dummyLong,lVec[mom],dummyLong2]+mass gamma[dummyLong,dummyLong2];
    gamma[dummyLong__,spinorV[{mom_,mass_,pol_}]]gamma[spinorVbar[{mom_,mass_,pol_}],dummyLong2__]^=gamma[dummyLong,lVec[mom],dummyLong2]-mass gamma[dummyLong,dummyLong2];
    spinorV[{mom_,mass_,pol_},sInd1_]spinorVbar[{mom_,mass_,pol_},sInd2_]^=gamma[sInd1,lVec[mom],sInd2]-mass delta[sInd1,sInd2] ;
	spinorU[{mom_,mass_,pol_},sInd1_]spinorUbar[{mom_,mass_,pol_},sInd2_]^=gamma[sInd1,lVec[mom],sInd2]+ mass delta[sInd1,sInd2];
	    res=Expand[expr];
	    res
]


simplifySpinChains[expr_,rules_]:=Block[{antiCom,sym,res,d},
(* anticommutation until canonical ordering is achieved *)
	gamma[sInd1_,s1VarLength___,ind1_,ind2_,s2VarLength___,sInd2_]/;Order[ind1,ind2]==-1:=-gamma[sInd1,s1VarLength,ind2,ind1,s2VarLength,sInd2]+antiCom[ind1,ind2] gamma[sInd1,s1VarLength,s2VarLength,sInd2];
	gamma[sInd1_,s1VarLength___,ind1_,ind2_,s2VarLength___,sInd2_]/;Order[ind1,ind2]==0:=gamma[sInd1,s1VarLength,s2VarLength,sInd2]sym[ind1,ind2];
	antiCom[sInd1_,lVec[sInd2_]]:=2 vector[sInd2,sInd1];
	antiCom[sInd1_,sInd2_]:=2 g[sInd2,sInd1];
	sym[lVec[sInd1_],lVec[sInd1_]]:=SP[sInd1,sInd1];
	sym[sInd1_,sInd1_]:=d;
	antiCom[lVec[sInd1_],lVec[sInd2_]]:=2 SP[sInd2,sInd1];
	gamma[sInd1_,sInd2_]:=delta[sInd1,sInd2]/;(!MemberQ[{spinorU,spinorV,spinorUbar,spinorVbar},Head@sInd1]&&!MemberQ[{spinorU,spinorV,spinorUbar,spinorVbar},Head@sInd2]);
res=expr //. rules;
res=contractLorentz[expr] //. rules
]


(* ::Input::Initialization:: *)
applyDiracEQ[expr_,rules_]:=Block[{antiCom,sym,res,permuteRight,permuteLeft,d},
(* anticommutation to the right until dirac eq. can be applied *)
permuteLeft[expr1_]:=Block[{gamma,spinorUbar,spinorVbar,res1,lVec},
gamma[spinorUbar[{mom_,mass_,pol_}],dummyLong1___,dummyShort_,lVec[mom_],dummyLong2__]:=-gamma[spinorUbar[{mom,mass,pol}],dummyLong1,lVec[mom],dummyShort,dummyLong2]+antiCom[dummyShort,lVec[mom]] gamma[spinorUbar[{mom,mass,pol}],dummyLong1,dummyLong2];
gamma[spinorUbar[{mom_,mass_,pol_}],lVec[mom_],dummyLong1__]:=mass gamma[spinorUbar[{mom,mass,pol}],dummyLong1];
gamma[spinorVbar[{mom_,mass_,pol_}],dummyLong1___,dummyShort_,lVec[mom_],dummyLong2__]:=-gamma[spinorVbar[{mom,mass,pol}],dummyLong1,lVec[mom],dummyShort,dummyLong2]+antiCom[dummyShort,lVec[mom]] gamma[spinorVbar[{mom,mass,pol}],dummyLong1,dummyLong2];
gamma[spinorVbar[{mom_,mass_,pol_}],mom_,dummyLong1__]:=-mass gamma[spinorVbar[{mom,mass,pol}],dummyLong1];
res1=expr1//Expand
    ];
permuteRight[expr1_]:=Block[{gamma,spinorU,spinorV,res1,lVec},
gamma[dummyLong1__,lVec[mom_],dummyShort_,dummyLong2___,spinorU[{mom_,mass_,pol_}]]:=-gamma[dummyLong1,dummyShort,lVec[mom],dummyLong2,spinorU[{mom,mass,pol}]]+antiCom[dummyShort,lVec[mom]] gamma[dummyLong1,dummyLong2,spinorU[{mom,mass,pol}]];
gamma[dummyLong1__,lVec[mom_],spinorU[{mom_,mass_,pol_}]]:=mass gamma[dummyLong1,spinorU[{mom,mass,pol}]];
gamma[dummyLong1__,lVec[mom_],dummyShort_,dummyLong2___,spinorV[{mom_,mass_,pol_}]]:=-gamma[dummyLong1,dummyShort,lVec[mom],dummyLong2,spinorV[{mom,mass,pol}]]+antiCom[dummyShort,lVec[mom]] gamma[dummyLong1,dummyLong2,spinorV[{mom,mass,pol}]];
gamma[dummyLong1__,lVec[mom_],spinorV[{mom_,mass_,pol_}]]:=-mass gamma[dummyLong1,spinorV[{mom,mass,pol}]];
res1=expr1//Expand
    ];

antiCom[sInd1_,lVec[sInd2_]]:=2 vector[sInd2,sInd1];
antiCom[sInd1_,sInd2_]:=2 g[sInd2,sInd1];
antiCom[lVec[sInd1_],lVec[sInd2_]]:=2 SP[sInd2,sInd1];
gamma[sInd1_,sInd2_]:=delta[sInd1,sInd2]/;(!MemberQ[{spinorU,spinorV,spinorUbar,spinorVbar},Head@sInd1]&&!MemberQ[{spinorU,spinorV,spinorUbar,spinorVbar},Head@sInd2]);

res=permuteRight[expr] //. rules;
    res=contractLorentz[res] //. rules;
res=permuteLeft[res] //. rules;
    res=contractLorentz[res] //. rules
]



performTraces[expr_,fileLoc_:$fileLocation]:=Block[{gammaTrace,tracesForm,tracesDef,g,lVec,myInd1,myInd2,ind1,ind2,vec1,vec2,SP,vector},
SetAttributes[g,Orderless];
SetAttributes[SP,Orderless];
tracesForm=FromFormToMathe[fileLoc];
tracesDef=Table[ToString[gammaTrace@@Table[ToExpression["traceInd"~~ToString[i]~~"_"],{i,j}]]~~"="~~ToString[tracesForm[[j]]],{j,Length@tracesForm}];
ToExpression/@tracesDef;
g[lVec[vec1_],lVec[vec2_]]=SP[vec1,vec2];
g[ind1_,lVec[vec1_]]:=Module[{ind},g[ind1,ind]vector[vec1,ind]];
expr
];


contractLorentz[expr_,opts:OptionsPattern[]]:=Block[{resultContracted,vector,g,d,delta,SP,gamma,gammaTrace,spinorU,spinorUbar,spinorV,spinorVbar,spinArg,vec1,vec2,vec3,ind1,ind2,ind3,s1Ind,s2Ind,s3Ind,s1VarLength,s2VarLength,aNum,lVec,tracesFile,addRules,simplSpinChain,spinSums},
    SetAttributes[SP,Orderless];
    SetAttributes[g,Orderless];
	SetAttributes[delta,Orderless];
	
	(* read oprtions *)
	addRules = OptionValue[additionalRules];
	tracesFile = StringJoin[OptionValue[traceFileLocation],"/trace_form.log"];
	simplSpinChain=OptionValue[spinChainSimplify];
	spinSums=OptionValue[polarizationSum];
	(* check input *)
	If[!FileExistsQ[tracesFile],
	Return["Error: "<>tracesFile<>" does not exist!"];
	Abort[]
	];
	
	If[!FreeQ[TrueQ[Head[#]==Rule||Head[#]==RuleDelayed]&/@addRules,False],
	Return["Error: additional rules contain non valid replacement!"];
	Abort[]
	];
	
	If[!BooleanQ[simplSpinChain],
	Print["Warning: Option spinChainSimplify is not a of type Boolean."];
	];
	
	If[!BooleanQ[spinSums],
	Print["Warning: Option polarizationSum is not a of type Boolean."];
	];
	
	(* Define algebra *)
	(* vectors *)
	vector[vec1_+vec2_,ind1_]=vector[vec1,ind1]+vector[vec2,ind1];
	vector[aNum_ vec1_,ind1_]/;NumericQ[aNum]=aNum vector[vec1,ind1];
	vector[vec1_,ind1_]vector[vec2_,ind1_]^=SP[vec1,vec2];
	Power[vector[vec1_,ind1_],aNum_]/;NumericQ[aNum]^=Power[SP[vec1,vec1],(aNum/2)];
(* scalar products *)
	SP[vec1_+vec2_,vec3_]=SP[vec1,vec3]+SP[vec2,vec3];
	SP[aNum_ vec1_,vec2_]/;NumericQ[aNum]=aNum SP[vec1,vec2];
(* metric *)
	g[ind1_,ind2_]*g[ind1_,ind3_]^=g[ind2,ind3];
	g[ind1_,ind1_]/;!NumericQ[ind1]=d;
	g[ind1_,ind2_]*vector[vec1_,ind2_]^=vector[vec1,ind1];
	Power[g[ind1_,ind2_],aNum_]/;NumericQ[aNum]^=Power[d,(aNum/2)];
	g[ind1_,ind2_]gamma[s1Ind__,ind1_,s2Ind__]^=gamma[s1Ind,ind2,s2Ind];
(*gamma *)
	gamma[s1VarLength__,ind1_,s2VarLength__]vector[vec1_,ind1_]^=gamma[s1VarLength,lVec[vec1],s2VarLength];
	gamma[s1VarLength__,s2Ind_]gamma[s2Ind_,s2VarLength__]^=gamma[s1VarLength,s2VarLength];
	delta[s1Ind_,s3Ind_] gamma[s1Ind_,s2VarLength__]^=gamma[s3Ind,s2VarLength];
    delta[s1Ind_,s3Ind_] gamma[s2VarLength__,s1Ind_]^=gamma[s2VarLength,s3Ind];
	gamma[s1Ind_,s2VarLength__,s1Ind_]=gammaTrace[s2VarLength];
(*gamma + spinor *)
	spinorUbar[spinArg_,s1Ind_] gamma[s1Ind_,s1VarLength__]^=gamma[spinorUbar[spinArg],s1VarLength];
	spinorVbar[spinArg_,s1Ind_] gamma[s1Ind_,s1VarLength__]^=gamma[spinorVbar[spinArg],s1VarLength];      
	spinorU[spinArg_,s1Ind_] gamma[s1VarLength__,s1Ind_]^=gamma[s1VarLength,spinorU[spinArg]];
	spinorV[spinArg_,s1Ind_] gamma[s1VarLength__,s1Ind_]^=gamma[s1VarLength,spinorV[spinArg]];
(* pure spinors *)
	spinorUbar[spinArg_,s1Ind_] delta[s1Ind_,s2Ind_]^=spinorUbar[spinArg,s2Ind];
	spinorVbar[spinArg_,s1Ind_] delta[s1Ind_,s2Ind_]^=spinorVbar[spinArg,s2Ind];
	spinorU[spinArg_,s1Ind_] delta[s1Ind_,s2Ind_]^=spinorU[spinArg,s2Ind];
	spinorV[spinArg_,s1Ind_] delta[s1Ind_,s2Ind_]^=spinorV[spinArg,s2Ind];
	
resultContracted=(expr//Expand)//.addRules;
resultContracted=performTraces[resultContracted,tracesFile]//Expand;
resultContracted=resultContracted//.addRules;

If[spinSums==True,
resultContracted=spinorSum[resultContracted];
resultContracted=performTraces[resultContracted,tracesFile]//Expand;
];

If[simplSpinChain==True,
	resultContracted=simplifySpinChains[resultContracted,addRules];
	resultContracted=applyDiracEQ[resultContracted,addRules];
];
resultContracted
];



contractColor[expr_,indFundamental_,indAdjoint_ ]:=Block[{delta,T,resultContracted,f,rulesForContraction,TF,Nc,reorderColor},
SetAttributes[delta,Orderless];
reorderColor[myexpr_, pats_List] :=
  Module[{h, rls},
    rls = MapIndexed[x : # :> h[#2, Replace[x, rls, -1]] &, pats];
    HoldForm @@ {myexpr /. rls} //. h[_, x_] :> x
  ];
rulesForContraction=Dispatch@{
	delta[a_,b_]f[c__]/;(!FreeQ[{c},a]&&!StringMatchQ[ToString@a,ToString@b]):>ReplaceAll[f[c],a->b],
	f[a_,b_,c_]:>(-I/TF Module[{n1,n2,n3},T[a,indFundamental[n1],indFundamental[n2]]T[b,indFundamental[n2],indFundamental[n3]]T[c,indFundamental[n3],indFundamental[n1]]-T[a,indFundamental[n1],indFundamental[n2]]T[c,indFundamental[n2],indFundamental[n3]]T[b,indFundamental[n3],indFundamental[n1]]]), (* TF=1/2 *)
	delta[a_,b_]T[c__]/;(!FreeQ[{c},a]&&!StringMatchQ[ToString@a,ToString@b]):>ReplaceAll[T[c],a->b],
	Times[T[a_,i_,j_],T[a_,k_,l_]]:>TF(delta[i,l]delta[j,k]-Power[Nc,-1]delta[i,j]delta[k,l]), (* fierz *)
	T[a_,i_,i_]:>0, (* traceless *)
	T[a_,l_,m_]T[b_,m_,l_]:>TF delta[a,b], (* single trace *)
	delta[a_,b_]delta[b_,c_]/;!StringMatchQ[ToString@a,ToString@c]:>delta[a,c],
	Power[delta[indAdjoint[x_],indAdjoint[y_]],2]/;!StringMatchQ[ToString@x,ToString@y]:>Power[Nc,2]-1,(* adjoint representation *)
	delta[indAdjoint[x_],indAdjoint[x_]]:>Power[Nc,2]-1, (* adjoint *)
	Power[delta[indFundamental[x_],indFundamental[y_]],2]/;!StringMatchQ[ToString@x,ToString@y]:>Nc, (* fundamental rep. *)
	delta[indFundamental[x_],indFundamental[x_]]:>Nc,(* fundamental *)
	Power[T[b_,c_,d_],2]:>TF(Power[Nc,2]-1),
	Power[f[b_,c_,d_],2]:>2TF(Power[Nc,2]-1)Nc

};

resultContracted=FixedPoint[ ReleaseHold@(ReplaceRepeated[reorderColor[Expand[#],{delta,T,f}],rulesForContraction])&,expr];
resultContracted //.Times[T[a_,i_,j_],T[b_,j_,l_]]/;(!StringMatchQ[ToString@a,ToString@b]&&!StringMatchQ[ToString@i,ToString@l]):>T[a,b,i,l]
];
