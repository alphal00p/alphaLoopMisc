(* ::Package:: *)

Print["The documented functions in this package are: \n ?getUVApproximant \n ?getUVprocess
 ----------------------------------------- 
 Needs the package LTDTools.m
"
];


getUVprocess::usage="get the process you need from qgraf.
	INPUT:
		- inParticles: list of incomming particle names (string) defined in the qgraph model 
			example syntax: {\"gluon\",\"gluon\"}
		- outParticles: list of outgoing particle names (string) defined in the qgraph model 
			example syntax: {\"higgs\",\"higgs\"}
		- loopNumber (Integer): number of loops
		- outputName (String): name of the qgraf output
		- propagatorCount (Association): assignement of the exact number of virtual propagators you want to consider. Anti-particle are propagators count as particle propagators
			example syntax: Association@{\"higgs\"\[Rule]0,\"gluon\"\[Rule]0,\"photon\"\[Rule]0,\"ghost\"\[Rule]0,\"eminus\"\[Rule]0,\"higgs\"\[Rule]0,\"u\"\[Rule]0,\"t\"\[Rule]4};
			Four tops/anti-tops are running in the loop
	Additional information:
		- the qgraf model consists currently of the following particles:
			{higgs, gluon, photon, ghost, ghostbar, eminus, eplus, u, ubar, t, tbar}
	OPTIONS:
		Options[getUVprocess]: {runCommand\[Rule]\"\",verbose\[Rule]True,FRulesLocation->\"./ct_graph_generation/minFeynRulesQEDQCD.m\"};
		runCommand is the system command to run qgraf. If \"\", then a suitable default bash command will be used.
		FRulesLocation: location of Feynman-rules
"
getUVapproximant::usage="get the analytic uv-ct
	INPUT: graphs or list of graphs as obtained from getUVprocess
	OPTIONAL INPUT: MaxEpsOrder\[Rule]0
"



qgraphGeneration[inParticles_List, outParticles_List,loopNumber_Integer,outputName_String,propagatorCount_Association,verbose_,command_]:=
	Block[{preamble,instates,outstates,loops,props,datfile},
		preamble="output='"<>outputName<>"';\nstyle='../../uvCT.sty';\nmodel='./minQCDQEDH.mod';\n";
		instates="in=";
		Do[
			instates=instates<>" "<>inParticles[[ii]]<>"[p"<>ToString@ii<>"],"
		,{ii,Length@inParticles}];
		instates=StringDrop[instates,-1]<>";\n";
		
		outstates="out=";
		Do[
			outstates=outstates<>" "<>outParticles[[ii]]<>"[p"<>ToString@(ii+Length@inParticles)<>"],"
		,{ii,Length@outParticles}];
		outstates=StringDrop[outstates,-1]<>";\n";

		loops="loops="<>ToString@loopNumber<>";\nloop_momentum=k;\noptions=onshell,notadpole,nosnail, onepi;\n";

		props="";
		Do[
			props=props<>"true=iprop["<>part<>","<>ToString@propagatorCount[[part]]<>","<>ToString@propagatorCount[[part]]<>"];\n"
		,{part,Keys@propagatorCount}];
		datfile=preamble<>instates<>outstates<>loops<>props;
		Export["./ct_graph_generation/qgraf.dat",datfile,"Text"];
		If[command===""&&verbose==True,
			Print["You will run the command:\n"<>"cd ./ct_graph_generation && rm -rf "<>outputName<>" a.out && gfortran qgraf-3.4.2.f && ./a.out"];
			Print["If the syntax does not work on your system, adjust the option 'runCommand'"]		
		];
		Run["cd ./ct_graph_generation && rm -rf "<>outputName<>" a.out && gfortran qgraf-3.4.2.f && ./a.out"];
		If[verbose==True,
			Print@datfile
		];
]


Protect[runCommand,verbose,FRulesLocation,particleContent];
Options[getUVprocess]={runCommand->"",verbose->True,FRulesLocation->"./ct_graph_generation/minFeynRulesQEDQCD.m",particleContent->ToString/@{higgs, gluon, photon, ghost, ghostbar, eminus, eplus, u, ubar, t, tbar}};
getUVprocess[inParticles_List, outParticles_List,loopNumber_Integer,outputName_String,propagatorCount_Association,opts:OptionsPattern[]]:=Block[{graphs,fullNumerator,lo},
(* check input *)
	Do[
		If[!MemberQ[OptionValue[particleContent],ii],
			Print["The in state: "<>ii<>" is not a particle in the model: "<>ToString@(Flatten@Keys@propagatorCount)];
			Abort[]
		]
	,{ii,inParticles}];
	Do[
		If[!MemberQ[OptionValue[particleContent],ii],
			Print["The out-state: "<>ii<>" is not a particle in the model: "<>ToString@(Flatten@Keys@propagatorCount)];
			Abort[]
		]
	,{ii,outParticles}];
	qgraphGeneration[inParticles, outParticles,loopNumber,outputName,propagatorCount,OptionValue[verbose],OptionValue[runCommand]];
	graphs=importGraphs["./ct_graph_generation/"<>outputName,sumIsoGraphs->False];
(*
Do[
lo=Association@("tree_level_numerator"\[Rule](v@@Cases[graphs[[gg,"numerator"]],Alternatives@@{_in,_out},Infinity])/. in[x___]\[RuleDelayed]x /. out[x___]\[RuleDelayed]x);
graphs[[gg]]=Append[graphs[[gg]],lo];
,{gg,Length@graphs}];
*)
	graphs=processNumerator[graphs,OptionValue[FRulesLocation],additionalRules->{scalarProp[x__]:>1}]
]



ClearAll[getUVApproximant]
Protect[MaxEpsOrder];
Options[getUVApproximant]={MaxEpsOrder->0};
getUVApproximant[graphs_,opts:OptionsPattern[]]:=Block[
	{sdd,mygraphs=Flatten@{graphs},propPos,loopBasis,expandPropagators,propCount,maxDegree,projector,projectedNumerator,translateToTopos,epsExpand,result},
	
	(* ------------------------------------------------------------- HELPER FUNCTIONS ------------------------------------------------------------*)
	(*wrapper for doing the uv tensor-projection at a given rank*)	
	projector[myrank_,loopMomentum_]:=Block[{projec,ind,coeff,tens,g,projectedTensor},
		If[OddQ[myrank],
			Return[0]
		];
		SetAttributes[g,Orderless];
		ind=Table[lMomInd[1][i],{i,myrank}];
		tens=(DeleteDuplicates@(Times@@@(Apply[g,#,{1}]&/@(ArrayReshape[#,{2,2}]&/@Permutations[ind]))));
		projec=translateToFeynCalc[(Inverse@(Contract[translateToFeynCalc[KroneckerProduct[tens,tens]]])).tens];
		translateToFeynCalc[tens.(Contract[translateToFeynCalc[projec*(Product[vector[loopMomentum,i],{i,ind}])]])]
	];
	(* plug in scalar topo and expand in epsilon*)
	epsExpand[num_,epsOrder_]:=Block[{topo1,eps,mUV,epsExpExpr},
		topo1[n_]:=(-1)^n I Pi^(4-2eps) Gamma[n+eps-2]/Gamma[n] 1/(mUV^2)^(n+eps-2);
		epsExpExpr=(ChangeDimension[num,4])/. D->4-2eps;
		epsExpExpr=Normal@(Series[epsExpExpr,{eps,0,epsOrder}]);
		If[epsExpExpr=!=0,
			epsExpExpr=Association@(Table[Association@(pow->(SeriesCoefficient[epsExpExpr,{eps,0,pow}]//FCE//Simplify)),{pow,Exponent[epsExpExpr,eps,Min],epsOrder}]),
			epsExpExpr=Association@(epsOrder->0);			
		];
		epsExpExpr
	
	];
	
	(* function for translating to scalar products to topologies: maybe rework for 2-loop needed(?) *)
	translateToTopos[uvgraph_,epsPower_]:=Block[{graphsUV=uvgraph[["uv_graphs"]],sps,props,dd,eqs,sols,mUV,tp,num},
		num=0;
		Do[
			sps=Cases[graph[["uv_topo"]],_Pair,Infinity]//Union;
			props=Array[dd,Length@graph[["uv_topo"]]];
			eqs=Thread@(Equal[props,(graph[["uv_topo"]]/. uvProp[a_,b_]:>a-mUV^2)]);
			tp=(ToExpression@("topo"<>ToString@(Length@graph[["uv_topo"]])));
			(* solve SPs in terms of denominators *)
			sols=Flatten@(Solve[eqs,sps]);
			num=num+(((graph[["numerator"]]//.sols)*tp@(graph[["prop_powers"]])//Expand)//. tp[x_List]dd[y_]:>tp@(ReplacePart[x, y -> x[[y]]-1])//.tp[x_List](dd[y_]^a_):>tp@(ReplacePart[x, y -> x[[y]]-a]));
			num=num/.tp[x_List]:>tp@@x;
			If[!FreeQ[num,Alternatives@@props],
				Print["Did not manage to express everything in terms of scalar topologies."];
				Print[num];
				Abort[]
			]
		,{graph,graphsUV}];
		num=epsExpand[num,epsPower];
		Association@("analytic_uv_ct"->num)
	];
	
	
	(* wrapper to expand numerators *)
	expandPropagators[uvgraph_]:=Block[{graph=uvgraph,uvexpandedNumerator,topo,uvAsso,props,uvProp,lambdaUV,uvExpansion,referenceTopo,uvG},
		uvAsso=graph[[{"uv_info"}]];
		graph=KeyDrop[graph,{"uv_info","analyticTensorCoeff"}];	
		graph=Append[graph,Association@("uv_graphs"->{})];
		Do[
			(* we dont take masses into account, since the CT will simply be the original diagram expanded around vanininshing momenta, and then we IR-regulate *)
			props=Cases[graph[["momentumMap"]],x_/;!FreeQ[x,Alternatives@@(Flatten@(Keys@expansion))]];
			(* uv-scaling *)
			referenceTopo=uvProp[SPD[#,#]//ExpandScalarProduct,1]&/@(Flatten@(Keys@expansion));
			props=(Times@@(uvProp[SPD[#,#]//ExpandScalarProduct,1]&/@props))/. Pair[Momentum[x_,D],Momentum[y_,D]]:>lambdaUV^Length@(DeleteCases[{x,y},Alternatives@@(Flatten@(Keys@expansion))]) Pair[Momentum[x,D],Momentum[y,D]] ;						
			(* expand for vanishing momenta *)
			uvExpansion=Normal@Series[props,{lambdaUV,0,(Values@expansion)[[1]]}];
			(* replace derivative of propagators *)
			uvExpansion=uvExpansion /.Derivative[a_,0][uvProp][x_,n_Integer]:>Factorial[a](-1)^a uvProp[x,n]^(a+1)/. lambdaUV->1;
			(*small sanity check of the number of UV-propagators*)
			If[Complement[Union@(Cases[uvExpansion,_uvProp,Infinity]),referenceTopo]=!={},
				Print["There is a mismatch between the reference topology and the UV-approximant"];
				Print["The reference topo is: "<>ToString@referenceTopo];
				Print["The expanded propagators are : "<>ToString@(Union@(Cases[uvExpansion,_uvProp,Infinity]))];
				Abort[]
				];			
			uvExpansion=CoefficientRules[uvExpansion,referenceTopo];
			(* overwrite uv-info *)
			Do[			
				graph[["uv_graphs"]]=Append[graph[["uv_graphs"]],
					Association@{Association@("loopMomenta"->Variables@(Flatten@(Keys@expansion))),
					Association@("uv_topo"->referenceTopo),
					Association@("prop_powers"->(Keys@uvExpansion)[[count]]),
					Association@("numerator"->(Values@uvExpansion)[[count]]*graph[["numerator"]]),
					Association@("momentumMap"->graph[["momentumMap"]])
				}]
			,{count,Length@uvExpansion}];
			(* project numerator *)
			Do[
				uvG=extractTensCoeff[graph[["uv_graphs"]][[uvCount]]];
				projectedNumerator=(Plus@@(DiracSimplify@(Contract[Times[#[[2]],projector[#[[1,1]],uvG[["loopMomenta",1]]]]&/@(uvG[["analyticTensorCoeff"]])])));
				graph[["uv_graphs"]][[uvCount]]=KeyDrop[graph[["uv_graphs"]][[uvCount]],{"momentumMap","numerator"}];
				graph[["uv_graphs"]][[uvCount]]=Append[graph[["uv_graphs"]][[uvCount]],Association@("numerator"->projectedNumerator)];
				
			,{uvCount,Length@graph[["uv_graphs"]]}]
		,{expansion,uvAsso}];
		graph
	];
	(* ------------------------------------------------------------- HELPER FUNCTIONS END ------------------------------------------------------------*)
	
	
	
	(* ------------------------------------------------------------- ACTUAL GENERATION ------------------------------------------------------------*)
	result=Table[
		If[!KeyExistsQ[graph,"analyticTensorCoeff"],
			graph=extractTensCoeff[graph]
		];
		loopBasis=Subsets[graph[["loopMomenta"]]]/. {}->Nothing;
		If[Length@loopBasis>1,
			Print["Only one-loop implemented so far"];
			Abort[]
		];
		
		graph=Append[graph,Association@("uv_info"-><||>)];
		maxDegree=Table[Max@(graph[["analyticTensorCoeff"]][[All,1,kk]]),{kk,Length@(graph[["loopMomenta"]])}];
		Do[
			propCount=Length@Position[graph[["momentumMap"]],x_/;!FreeQ[x,Alternatives@@ll],1];
			sdd=(4*Length@ll-2*propCount + Total@@(Extract[maxDegree,Position[graph[["loopMomenta"]],x_/;!FreeQ[x,Alternatives@@ll],Heads->False]]));			
			graph[["uv_info"]]=Append[graph[["uv_info"]],Association@(ll->sdd)];				
		,{ll,loopBasis}];
		graph=expandPropagators[graph];
		graph=translateToTopos[graph,OptionValue[MaxEpsOrder]]
	,{graph,mygraphs}];
	If[Length@result==1,
		result[[1]],
		result	
	]
]
