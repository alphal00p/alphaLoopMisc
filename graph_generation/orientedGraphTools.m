(* ::Package:: *)

Print["The documented functions in this package are: \n ?plotGraph \n ?findIsomorphicGraphs \n ?constructCuts \n ?importGraphs \n ?getLoopLines \n ?getCutStructure \n ?writeMinimalJSON \n ?contractLorentz \n ?contractColor \n ?extractTensCoeff \n ?getSymCoeff \n ?processNumerator
 ----------------------------------------- 
 Needs the package IGraphM which can be downloaded from https://github.com/szhorvat/IGraphM. !!! \n Run: Get[\"https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m\"] for installation.
 Make sure you have the \"trace_form.log\" file for the gamma-algebra \n "
];
Needs["IGraphM`"];
constructCuts::usage="Finds all cutkovsky cuts of a graph.
	Input: 
		graph (see output from qgraf \"orientedGraph.sty\")
	Optional Input:
		NumberFinalStates (default:All): number of cuts
		DisplayCuts (default: True): plots graphs with high-lighted cuts
		DetectSelfEnergy (default: True): append entry to graph, which contains the information on the self energy insertions
	Output: all possible cut graphs
"
findIsomorphicGraphs::usage="Finds isomorphisms between graphs.
	Input: 
		List of graphs (see output from qgraf \"orientedGraph.sty\")
	Optional Input:
		edgeTagging (default: {}): List of Keyword(s) of association(s) in the graphs on which to tag by (e.g. \"partilcleType\")
	Output: List where every entry has the form
	{
		{Position of graph to be mapped, Positon of graph it maps to},
		{
		{isomorphism, permutation of edges, orientation changes: +1 = no orientation change, -1: orientation change}
		}
	}
"
plotGraph::usage="plots graph
		Input: 
			graph (see output from qgraf \"orientedGraph.sty\")
		Optional Input:
			edgeLabels(default: {}): List of Keyword(s) of association(s) in the graphs on which to tag by (e.g. \"partilcleType\")
			plotSize (default: Scaled[0.5])
"
importGraphs::usage="Import QGraf graphs generated with orientedGraphs.sty
			Input: 
				\"PathToQGrafOutput\"
			Optional Input:
				sumIsoGraphs\[Rule]True/False (* Default: True *): Align momenta of isomorphic scalar topologies and sum numerators
"
getLoopLines::usage="Appends loop-lines to graph.
		Input: 
			graphs (List or single graphs) (see: importGraphs[output from qgraf \"orientedGraph.sty\"])
		Output:
			graphs with appended association for loop-lines
"
getCutStructure::usage="Appends cut-strucuture to graph.
		Input: 
			graphs (List or single graphs) (see: importGraphs[output from qgraf \"orientedGraph.sty\"])
				- if graphs have no attribute \"loopLines\", it will be generated
				- if graphs have attribute \"contourClosure\" ->{\"Above\",...,\"Below\",..} this particular closure will be used. Otherwise
				  default closure from above.
		Output:
			graphs with appended association for cut-structure
"
writeMinimalJSON::usage="Exports JSON-File for LTD.
		Input: 
			graphs (List or single graphs) (see: importGraphs[output from qgraf \"orientedGraph.sty\"])
				- if graphs have no attribute \"cutStructure\", it will be generated
				- if graphs have no attribute \"name\", name will be genrated as PROCESSNAME_DIAG_DIAGNUMBER
			numAssociation (type: association): association of masses and external momenta to numerical values
		Optional Input:
			Options: {processName\[Rule] \"NameOfProcess\"(Default: \"\"), exportDirectory\[Rule] \"PathToExportDir/\" (Default: \"./\"), exportDiagramwise (Default: False),writeNumerator (Default: False)}
		Output:
			None"
contractLorentz::usage="Contracts Lorentz-indices and Spinor indices and performs Gamma-traces.: 
	 It works with the following objects:
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
"
contractColor::usage=" Performs the color-algebra. It works with the following objects:
		delta[sunA[x1],sunA[x2]],delta[sunF[x1],sunF[x2]] : Kronecker-delta k,n either \"adjoint-\" or \"fundamental\"-indices
		f[sunA[x1],sunA[x2],sunA[x3]]: structure constant of su(Nc) (sunA[x_]): \"adjoint-indices\"
		T[sunA[a],sunF[l],sunF[m]]: (l,m)-th entry of the a-th generator of fundamental representation with (a): \"adjoint-index \" and components (l,m): \"fundamental-indices\"
	INPUT: 
		expr: The expression on which color algebra should be performed		
	OUTPUT: Expression with performed color-algebra.
	NEW CONSTANTS: Nc from su(Nc), TF from trace normalization of su(Nc) (e.g. TF=1/2)
"
extractTensCoeff::usage=" Extracts from a given numerator the loop-momentum tensor coefficients and returns them ordered by rank.
	INPUT: 
		graph (or list of graphs): see output from importGraphs[], needs attribute  \"numerator\" .
	OPTIONAL INPUT: 
		consistencyCheckLevel (Default=1): INTEGER: 2: Run inbetween checks (for debugging);
											  1: Perform Dirac-/Lorentz-algebra and compare final result vs input expression
											  0: No-checks
	OUTPUT: graph with appended entry \"analyticTensorCoeff\" (see example)
	Example:
	graph=<|...,\"numerator\"->SP[k1,p1]SP[p1,p2]+SP[p1,p2]+SP[k1,k1]+SP[k1,k2] |>(* numerator , 2-loop *);
	extractTensCoeff[graph]				
	OUTPUT (example: ):
	<|
     ...
	\"analyticTensorCoeff\"\[Rule]{
      {{0,0},SP[p1,p2]}, (* rank 0 *)
      {{1,0},SP[p1,p2] vector[p1,lMomInd[1][1]]}, (* rank one in k1, lMomInd[X][Y]]=k_X^{mu[Y]} *)
      {{2,0},g[lMomInd[1][1],lMomInd[1][2]]},{{1,1},g[lMomInd[1][1],lMomInd[2][1]]} (*{{k1^mu1 k1^mu2},...},{{k1^mu1 k2^nu1},...}*)
     }
    |>
"
getSymCoeff::usage=" Gives the symmetric sum over tensor coefficients with explicitly expanded numeric Lorentz indices. 
	INPUT: 
		graph (or list of graphs) with attribute \"numerator\" and/or \"analyticTensorCoeff\". 
    OPTIONAL INPUT
        outputFormat (default: \"long\"): 
             - \"short\": only non-vanishing coefficients are produced with appended index information.
             - \"long\":  all coefficients are produced in canonical order. No index information is appended.
	OUTPUT: graph with appended entry \"symmetrizedExpandedTensorCoeff\"
"
processNumerator::usage="Processes the numerator of a graph.
  INPUT:
   - graph (or list of graphs)
   - absolute path to Feynamn rules
  OPTIONAL INPUT:
   - algebraOptions\[Rule]{additionalRules -> {1->1}, traceFileLocation->$fileLocation,spinChainSimplify->False,polarizationSum->False } : see ?contractLorentz
   - extractTensors\[Rule]False/True (True: analytic tensor coefficients are appended : see ?extractTensCoeff )
   - symCoefficients\[Rule]True/False ( True:  symmetric sum over tensor coefficients with explicitly expanded numeric Lorentz indices are appendend:  ?getSymCoeff )
   - coeffFormat\[Rule]\"long\"/\"short\" : Specifier of symmetric coefficient format
  OUTPUT:  Graph with processed numerator
"



Options[importGraphs]={sumIsoGraphs->True};
importGraphs[file_:String,opts:OptionsPattern[]]:=Block[{graphs, loopMom,isos,giso,kNew,relabelLoopMom,graphToBeMapped,graphToMapTo,momToMapTo,momToBeMapped,solMom},

graphs=ToExpression["{"<>StringDrop[StringReplace[Import[file,"Text"],{"**NewDiagram**"->",","\n"->""}]<>"}",1]];
loopMom=Complement[Variables@(Union@Flatten@graphs[[All,"momentumMap"]]/.f_[x_]:>x),Union@(Cases[graphs[[All,"momentumMap"]],f_[x_]:>x,Infinity] )]
;
graphs=Append[#,"loopMomenta"->loopMom]&/@graphs;
(* only relevant for forward scattering *)
graphs=DeleteCases[graphs,x_/;!FreeQ[x[["momentumMap"]],0],{1}];
If[OptionValue[sumIsoGraphs]==True,
isos=findIsomorphicGraphs[graphs,edgeTagging->{"particleType"}];
giso=Length@isos;
(* filter isomorphisms *)

Do[
isos=Join[isos[[1;;countIso]],DeleteCases[isos[[countIso+1;;]],x_/;!FreeQ[x[[1]],isos[[countIso,1,1]] ],{1}]];
If[countIso==Length@isos-1,
Break[]
]
,{countIso,giso}];
(* build minimal graphs *)
relabelLoopMom=Array[kNew,Length@loopMom];
Do[

graphToBeMapped=graphs[[iso[[1,1]]]];
graphToMapTo=graphs[[iso[[1,2]]]];
momToMapTo=graphToMapTo[["momentumMap"]];
momToBeMapped=Permute[graphToBeMapped[["momentumMap"]],First[Sort[iso[[2,All,"permutation"]]]]]/. Thread[Rule[loopMom,relabelLoopMom]];
(* we consider the isomorphism as a pure shift in the loop-momenta which mappes the scalar subtopologies *)
solMom=(Solve[(momToMapTo^2-momToBeMapped^2)==0,relabelLoopMom]);
If[Length@solMom==0,
Print["Error in alignment of momenta. Please run with: ``sumIsoGraphs\[Rule]False``"];
Abort[],
solMom=solMom[[1]]
];
graphs[[iso[[1,2]],"numerator"]]=graphs[[iso[[1,2]],"numerator"]]+(graphToBeMapped[["numerator"]]//. Thread[Rule[loopMom,relabelLoopMom]]//. solMom);
graphs[[iso[[1,1]]]]=Append[graphs[[iso[[1,1]]]],"delete"->{True,iso[[1]]}];
,{iso,isos}];
(* delete mapped graphs *)
graphs=DeleteCases[graphs,x_/;KeyExistsQ[x,"delete"],{1}];
(* last check for aligned scalar topologies: check that every summand in every numerator has the same scalar props.... maybe a bit overcautious? *)
Do[
If[Length@((Sort/@(Cases[#,_prop,Infinity]&/@(Flatten@Apply[List,(Flatten@{gg[["numerator"]]}),{1}])/. prop[f_[x_,y_]]:>x^2//FullSimplify))//Union)>1,
Print["Error in alignment of momenta. Please run with: ``sumIsoGraphs\[Rule]False``"];
Abort[]
]
,{gg,graphs}]
];
graphs
]



Protect[NumberFinalStates,DisplayCuts,DetectSelfEnergy];
Options[constructCuts] ={ NumberFinalStates->All,DisplayCuts->True,DetectSelfEnergy-> True};
constructCuts[graph_,opts:OptionsPattern[]]:=Block[{finalStateNum,detectSE,displayCuts,particles,momenta, trees,edges,cutEdges,forests,external,vertices,tmpDiag,cutTag,loopNum,cutGraphs,seInfo,props,cutProps,cutPos,prop,seEdges,probEdge,selfEnergy,subGEdges,sePostion,sePos,findAllSpanningTrees,externalMom},
SetAttributes[UndirectedEdge,Orderless];
(* explorative function to construct all spanning trees *)
findAllSpanningTrees[edges_,external_,loopNumber_]:=
	Block[{internalEdges,spanningTree,partitions,openEdges,addEdge},
	(* add edge connected to tree *)
		addEdge[openEdge_,tree_]:=Block[{shortEdges,addTree},
				shortEdges=List/@Cases[openEdge,UndirectedEdge[a_,b_]/;(MemberQ[VertexList@(Graph@tree),a]&&!MemberQ[VertexList@(Graph@tree),b])||(!MemberQ[VertexList@(Graph@tree),a]&&MemberQ[VertexList@(Graph@tree),b]),Infinity];
				If[Length@shortEdges>0,
					addTree=Join[tree,#]&/@shortEdges,
					addTree={tree}]
				];
		internalEdges=Complement[edges,external];
		spanningTree={{internalEdges[[1]]}};
		(* add edges until all trees are found *)
		Do[
			openEdges=Complement[internalEdges,#]&/@spanningTree;
			spanningTree=Flatten[Table[addEdge[openEdges[[treeCount]],spanningTree[[treeCount]]],{treeCount,Length@openEdges}],1];
		,{kk,Length@internalEdges-loopNumber}];
		Join[external,#]&/@((Union@(Sort/@spanningTree)))
	];

(* detects selfenery insertions *)
selfEnergy[edgesSubG_,problematicEdge_,extEdges_]:=
	Block[{graphs,result},
		graphs=Complement[edgesSubG,problematicEdge];(* two-component graph*)
		graphs=EdgeList/@(WeaklyConnectedGraphComponents@graphs);
		(* the component not having externals is the SE *)
		If[ContainsAny[graphs[[1]],extEdges],
			result=graphs[[2]],
			result=graphs[[1]]
		];
		result
	];

(* initialize data *)
If[!MemberQ[Keys@graph,"momentumMap"],
	Print[" \"momentumMap\" missing"];
	Abort[]
];
If[!MemberQ[Keys@graph,"edges"],
	Print[" \"edges\" missing"];
	Abort[]
];
If[!MemberQ[Keys@graph,"particleType"],
	Print[" \"particleType\" missing"];
	Abort[]
];
If[!MemberQ[Keys@graph,"loopMomenta"],
	Print[" \"loopMomenta\" missing. Import with function importGraphs[]"];
	Abort[]
];

finalStateNum=OptionValue[NumberFinalStates];
displayCuts=OptionValue[DisplayCuts];
detectSE=OptionValue[DetectSelfEnergy];

edges=graph[["edges"]]/. Rule[a_,b_]:>UndirectedEdge[a,b] /. DirectedEdge[a_,b_]:>UndirectedEdge[a,b]/. TwoWayRule[a_,b_]:>UndirectedEdge[a,b];
momenta=graph[["momentumMap"]];
particles=graph[["particleType"]];

externalMom=Cases[momenta,x_/;FreeQ[x,Alternatives@@(Join[graph[["loopMomenta"]],graph[["loopMomenta"]]])]]//Union;

external=Extract[edges,Position[graph[["momentumMap"]],Alternatives@@externalMom]];

vertices=VertexList@edges;
loopNum=Length@(graph[["loopMomenta"]]);
If[Length@(WeaklyConnectedGraphComponents@Complement[edges,external])>1,
Return[{Append[Append[graph,"cutInfo"->"factorizedLoopDiagram"],"selfEnergieInfo"->"factorizedLoopDiagram"]}];
Break[]
];

trees=findAllSpanningTrees[edges,external,loopNum];
(* deleting edge from spanning tree gives forest *)
forests=DeleteDuplicates@(Sort/@(Join[external,#]&/@(Flatten[Subsets[Complement[#,external],{Length@#-1-Length@external}]&/@(trees),1])));

(* a two-forest must have all vertices of the original graph: if it has not, the cut is improper \[Rule] delete all of them *)
forests=DeleteCases[forests,x_/;Length@(VertexList@x)!=Length@vertices,{1}];

(* all edges not in the 2-forest are cut *)
cutEdges={Flatten@(EdgeList/@WeaklyConnectedGraphComponents[#,in[x_]]),Flatten@(EdgeList/@WeaklyConnectedGraphComponents[#,out[x_]]), Complement[edges,#]}&/@forests;

(* Filtering *)
(* 1.) A component having ''in'' as well as ''out'' states is not allowed! *)
cutEdges=DeleteCases[cutEdges,x_/;!FreeQ[x[[1]],out[a_]],{1}];
cutEdges=DeleteCases[cutEdges,x_/;!FreeQ[x[[2]],in[a_]],{1}];

(* 2.) A cut edge is not allowed to be complete in the right or left diagram. If it is, uncut it *)
cutEdges=
	Table[
		tmpDiag=cutDiag;
		Do[
			If[(ContainsAll[VertexList@cutDiag[[1]],VertexList@({edge})]),
			tmpDiag[[1]]=Join[tmpDiag[[1]],{edge}];
			tmpDiag[[-1]]=Complement[tmpDiag[[-1]],{edge}];
		];
		If[(ContainsAll[VertexList@cutDiag[[2]],VertexList@({edge})]),
			tmpDiag[[2]]=Join[tmpDiag[[2]],{edge}];
			tmpDiag[[-1]]=Complement[tmpDiag[[-1]],{edge}];
		];
		,{edge,cutDiag[[-1]]}];
		Sort/@tmpDiag
	,{cutDiag,cutEdges}];
cutEdges=DeleteDuplicates@cutEdges;

(* Add assuciation cut-graph *)
cutTag=Table[edges /. Thread@Rule[cutEdges[[cutCount,1]],"left"]/.Thread@Rule[cutEdges[[cutCount,2]],"right"] /. Thread@Rule[cutEdges[[cutCount,3]],"cut"],{cutCount,Length@cutEdges}];

If[NumberQ[finalStateNum],
	cutTag=DeleteCases[cutTag,x_/;Count[x,"cut"]!=finalStateNum,{1}]
];

If[displayCuts==True,
	Do[
		Print@HighlightGraph[edges,Extract[edges,Position[cutTag[[cc]],"cut",Heads->False]]]
	,{cc,Length@cutTag}]
];
cutGraphs=Table[Append[graph,"cutInfo"->cutTag[[ct]]],{ct,Length@cutTag}];

If[detectSE==True,
	seInfo=Table[
		cutPos=Position[gg[["cutInfo"]],"cut"];
		props=prop[#^2//Expand //FullSimplify]&/@gg[["momentumMap"]];
		cutProps=Extract[props,cutPos];
		(* info on edges which will blow up *)
		probEdge=(Extract[#,Complement[Position[props,Alternatives@@cutProps],cutPos]]&/@{gg[["edges"]],gg[["particleType"]],gg[["momentumMap"]],gg[["cutInfo"]]});
		(* get the momenta of the selfenergy insertions *)
		seEdges=Table[
			subGEdges=Extract[gg[["edges"]],Position[gg[["cutInfo"]],probEdge[[-1,countInsertions]]]];
			selfEnergy[subGEdges,{probEdge[[1,countInsertions]]},external]
		,{countInsertions,Length@probEdge[[-1]]}];
		If[Length@(Flatten@probEdge)==0,
			probEdge={}
		];
		<|"selfEnergieInfo"-><|"problematicEdges"->probEdge ,
		"seDiagrams"->Table[
			sePos=Position[edges,Alternatives@@se];
			<|"EdgesSE"->Extract[graph[["edges"]],sePos],
			"particleTypeSE"->Extract[gg[["particleType"]],sePos],
			"momentumMapSE"->Extract[gg[["momentumMap"]],sePos],
			"cutInfoSE"-> Extract[gg[["cutInfo"]],sePos]
		|>
	,{se,seEdges}]|>|>
	,{gg,(cutGraphs/. Rule[a_,b_]:>UndirectedEdge[a,b] /. DirectedEdge[a_,b_]:>UndirectedEdge[a,b]/. TwoWayRule[a_,b_]:>UndirectedEdge[a,b])}];
	cutGraphs=Table[Append[cutGraphs[[gg]],seInfo[[gg]]],{gg,Length@cutGraphs}]
];
	cutGraphs
]


(* represent graph by: {{edge1: (Rule[a,b] for directed edge, TwoWayRule[a,b] for undirected edge),  ,tag1:List},...,{edgeN,tagN:List}} *)
Protect[edgeTagging];
Options[findIsomorphicGraphs] ={ edgeTagging->{}};
findIsomorphicGraphs[graphs_:List,tag:OptionsPattern[]]:=Block[{edges,myTaggedGraphs,tagKeys,tagging,taggingToNum,$uniqueTag, graphsEdgeColored,graphsEdgeColoredSimple,isIsomorph,$myAssoc,isomorphisms,shifts,g1,g2,shiftsAndPermutations,orientation,mappedGraph,graphID1,graphID2,permutation,mapping},
	SetAttributes[TwoWayRule,Orderless];
	edges=graphs[[All,"edges"]];
	If[!FreeQ[edges,_Missing],
		Print["graph is missing Key :\"edges\""];
		Abort[]		
		];
	tagKeys = OptionValue[edgeTagging];

	(* taggings become colors on edge-colored graphs, multiple taggings will be translated to single prime which becomes the color of the edge *)
	If[Length@tagKeys!=0,
		tagging=Values/@graphs[[All,tagKeys]];
		If[Length@(Union@(Length/@tagging))!=1||Length@(Union@Map[Length,tagging,{2}])!=1,
			Print["Not every edge has a tagging for choosen edgeTagging: "<>ToString[tagKeys]];
			Abort[]		
		];
		If[!FreeQ[tagging,_Missing],
			Print["Not every edge has a tagging for choosen edgeTagging: "<>ToString[tagKeys]];
			Abort[]		
		];
		tagging=Transpose/@tagging,		
		(* if there is not tagging, introduce random tagging *)
		tagging=(edges/. Rule[a_,b_]:>{$uniqueTag} /. TwoWayRule[a_,b_]:> {$uniqueTag}/. DirectedEdge[a_,b_]:> {$uniqueTag}  /. UndirectedEdge[a_,b_]:> {$uniqueTag});
		];
	
(* start prime range late, so that translation of multi-edge graphs is still unique *)
	taggingToNum=Thread@(Rule[Variables@tagging,Prime[Range[100,Length@(Variables@tagging)+99]]]);

(* translate to undirected graphs and apply tagging *)
	myTaggedGraphs=MapThread[Riffle[{#1},{#2}]&,{edges,tagging},2];

	graphsEdgeColored=myTaggedGraphs /. taggingToNum /. {DirectedEdge[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x} /. {UndirectedEdge[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x}/. {Rule[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x} /. {TwoWayRule[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x};

(* get isomorphism for graphs which are isomorphic; IGVF2 needed for multi-edge graphs: 
	There is no algorithm working for multi-edge graphs directly. Therefore translate into edge-colored simple graphs! The color is the multiplicity of the edge e.g.: g1Multiedge={{1\[TwoWayRule]4,5},1\[TwoWayRule]4,5}} --> g1ColouredSimple={1\[TwoWayRule]4\[Rule]2*5}*)
	graphsEdgeColoredSimple=(Counts/@(Sort/@graphsEdgeColored)) /.Association:>$myAssoc/.Rule[List[TwoWayRule[a_,b_],color1_],color2_]:> Rule[TwoWayRule[a,b],color1*color2] /. TwoWayRule->TwoWayRule /. $myAssoc:>Association;

(* thats a bit of laziness, since over computation of half the entries.... *)
	isIsomorph=Table[IGVF2IsomorphicQ[{Graph[VertexList[Keys@graph1],Keys[graph1]],"EdgeColors"->graph1},{Graph[VertexList[Keys@graph2],Keys[graph2]],"EdgeColors"->graph2}],{graph1,graphsEdgeColoredSimple},{graph2,graphsEdgeColoredSimple}]//LowerTriangularize;
	isomorphisms=(Position[isIsomorph,True]/. {i_,i_}-> Nothing);

(* compute isomorphism for isomorphic graphs *)
	shifts=ConstantArray[0,Length[isomorphisms]];
	Do[
		g1=graphsEdgeColoredSimple[[isomorphisms[[i,1]]]];
		g2=graphsEdgeColoredSimple[[isomorphisms[[i,2]]]];
		shifts[[i]]={{isomorphisms[[i,1]],isomorphisms[[i,2]]},IGVF2FindIsomorphisms[{Graph[VertexList[Keys@g1],Keys[g1]],"EdgeColors"->g1},{Graph[VertexList[Keys@g2],Keys[g2]],"EdgeColors"->g2}]};
	,{i,1,Length[isomorphisms]}];
(* find edge permutations *)
	shiftsAndPermutations=Table[
		Join[{shifts[[countShifts,1]]},
		{Table[
		
			{shifts[[countShifts,2,isoCount]],
			If[Length@(Complement[graphsEdgeColored[[shifts[[countShifts,1,1]]]]/.shifts[[countShifts,2,isoCount]],graphsEdgeColored[[shifts[[countShifts,1,2]]]]])==0,
			FindPermutation[graphsEdgeColored[[shifts[[countShifts,1,1]]]]/.shifts[[countShifts,2,isoCount]],graphsEdgeColored[[shifts[[countShifts,1,2]]]]],
			"errorPerm"
			]}
		,{isoCount,Length@shifts[[countShifts,2]]}]}]
	,{countShifts,Length@shifts}];
	
	shiftsAndPermutations=shiftsAndPermutations /. {x_,"errorPerm"}:>Nothing;
	
(* find if edge direction changed: No=1,Yes=-1*)
	Do[
		graphID1=shiftsAndPermutations[[shiftCount,1,1]];
		graphID2=shiftsAndPermutations[[shiftCount,1,2]];
		Do[
			mapping=shiftsAndPermutations[[shiftCount,2,isoCount,1]];
			permutation=shiftsAndPermutations[[shiftCount,2,isoCount,-1]];
			mappedGraph=Permute[myTaggedGraphs[[graphID1]],permutation]/.mapping /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b});
			orientation=(Keys[mappedGraph[[All,1]]]-Values[mappedGraph[[All,1]]])/(Keys[(myTaggedGraphs[[graphID2,All,1]] /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b}))]-Values[(myTaggedGraphs[[graphID2,All,1]] /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b}))])//Simplify;
			shiftsAndPermutations[[shiftCount,2,isoCount]]=Join[shiftsAndPermutations[[shiftCount,2,isoCount]],{orientation}];
		,{isoCount,Length@shiftsAndPermutations[[shiftCount,2]]}];
	,{shiftCount,Length@shiftsAndPermutations}];
	shiftsAndPermutations=Table[{shiftsAndPermutations[[countIso,1]],Association/@Table[Thread@Rule[{"vertexMap","permutation","orientationChange"},iso],{iso,shiftsAndPermutations[[countIso,2]]}] },{countIso,Length@shiftsAndPermutations}]
];


(* GraphComputation`GraphPlotLegacy is a hack, since mathematica updated GraphPlot *)
Protect[edgeLabels,plotSize]
Options[plotGraph] ={ edgeLabels->{},plotSize->Scaled[0.5]};
plotGraph[graph_,opt:OptionsPattern[]]:=Block[{mygraph,imSize=OptionValue[plotSize]},
	Do[
	mygraph=Transpose@(Values@(gg[[Join[{"edges"},Flatten@{OptionValue[edgeLabels]}]]]))/. TwoWayRule->Rule /. UndirectedEdge->Rule /. DirectedEdge->Rule /. {Rule[a_,b_],c__}:>{Rule[a,b],Flatten@{c}} /. {Rule[a_,b_]}:>{Rule[a,b],{Rule[a,b]}};
	Print[GraphComputation`GraphPlotLegacy[mygraph,EdgeRenderingFunction->({{RGBColor[0.22,0.34,0.63],Arrowheads[0.015],Arrow[#]},If[#3=!=None,Text[ToString[#3],Mean[#1],Background->White],{}]}&),MultiedgeStyle->.3,SelfLoopStyle->All,ImageSize->imSize]];
	,{gg,Flatten@{graph}}]
];



getLoopLines[graphs_]:=Block[{eMom,lMom,lLines,lProp,orientation,res,mass,pos},
res=Table[
If[!KeyExistsQ[graph,"loopMomenta"],
	Print["Error: Graphs have no entry \"loopMomenta\". Import with importGraph[]"];
	Abort[]
];
eMom=Complement[Variables@(Union@Flatten@graph[["momentumMap"]]),graph[["loopMomenta"]]];
lLines=CoefficientArrays[((graph[["momentumMap"]] /. Thread[Rule[eMom,0]]/. 0:>Nothing)^2 //FullSimplify)/.(x_)^2:>x// DeleteDuplicates,graph[["loopMomenta"]]][[2]]//Normal;

lProp=Table[
	orientation=(ReplaceAll[#,x_/;!NumberQ[x]:>0]&/@(FullSimplify@((graph[["momentumMap"]]/. Thread[Rule[eMom,0]])/(Dot[graph[["loopMomenta"]],ll])) ));
	pos=Position[Abs@orientation,1,Heads->False];
	<|"signature"->ll,
	"propagators"->Table[
	<|"mass"->mass[graph[["particleType",pp]]],
	  "name"-> If[KeyExistsQ[graph,"name"],
	  graph[["name",pp]],
	  "prop"<>ToString@pp
	  ],
	  "power"->If[KeyExistsQ[graph,"power"],
	  graph[["power",pp]],
	  1
	  ],
	  "q"->(graph[["momentumMap",pp]]*orientation[[pp]] /. Thread[Rule[graph[["loopMomenta"]],0]] /. 0->ConstantArray[0,4])
	  |>
	  ,{pp,Flatten@pos}]	  
	|>	 
,{ll,lLines}];
Append[graph,<|"loopLines"->lProp|>]
,{graph,Flatten@{graphs}}];
If[Length@res==1,
res=res[[1]]
];
res
]




$fileDirectory=If[$Input==="",FileNameJoin[{NotebookDirectory[],""}],FileNameJoin[{DirectoryName@$InputFileName,""}]];
getCutStructure[graphs_]:=Block[{session,result,pyTest,myGraphs,ctStruc,ll},

If[!FileExistsQ[$fileDirectory<>"/compute_cut_structure.py"],
	Print["No cut-structure script found in: \""<>$fileDirectory<>"/compute_cut_structure.py"<>"\""];
	Abort[]
];
session=StartExternalSession[<|"System"->"Python","Version"->"3","ReturnType"->"String"|>];
If[!StringMatchQ[session["Version"],"3"~~___],
	Print["Error: No external python3 evaluator registered."];
	Print["see: \" https://reference.wolfram.com/language/workflow/ConfigurePythonForExternalEvaluate.html \""];
	Abort[]
];
pyTest=ExternalEvaluate[session,"import sys"];

If[pyTest=!="Null",
	Print["Error in ExternalEvaluate[session,\"import sys\"]"];
	Print["possible cause: Usage of python 3.8"];
	Print["possible fix: \n Change \" yourPathTo/WolframClientForPython/wolframclient/utils/externalevaluate.py \""];
	Print["66c66
		<         exec(compile(ast.Module(expressions, []), '', 'exec'), current)
		---
		>         exec(compile(ast.Module(expressions), '', 'exec'), current)
	"];
	Abort[]
];
ExternalEvaluate[session,"sys.path.append('"<>$fileDirectory<>"')"];
ExternalEvaluate[session,"import compute_cut_structure as cct"];

If[ContainsAny[KeyExistsQ[#,"loopLines"]&/@(Flatten@{graphs}),{False}],
	myGraphs=getLoopLines[graphs],
	myGraphs=graphs
];
result=Table[
	ll=gg[["loopLines"]][[All,"signature"]];
	ExternalEvaluate[session,"loop_line_signatures = "<>ExportString[ll,"PythonExpression"]];
	ExternalEvaluate[session,"cut_structure_generator = cct.CutStructureGenerator(loop_line_signatures)"];
	If[KeyExistsQ[gg,"contourClosure"],
		ExternalEvaluate[session,"contour_closure = "<>StringReplace[ExportString[(gg[["contourClosure"]]/. "Above"->"cct.CLOSE_ABOVE"/. "Below"->"cct.CLOSE_BELOW"),"PythonExpression"],"\""->"" ]],
		ExternalEvaluate[session,"contour_closure = "<>StringReplace[ExportString[ConstantArray["cct.CLOSE_ABOVE",Length@ll[[1]]],"PythonExpression"],"\""->"" ]];
	];
	ctStruc=ImportString[ExternalEvaluate[session,"cut_structure_generator(contour_closure)"],"PythonExpression"];
	Append[gg,"cutStructure"->ctStruc]
,{gg,Flatten@{myGraphs}}];
Clear[session];
If[Length@result==1,
	result=result[[1]]
];
result
]



Protect[processName,exportDirectory,exportDiagramwise,writeNumerator];
Options[writeMinimalJSON]={processName->"",exportDirectory->"./",exportDiagramwise->False,writeNumerator->False}
writeMinimalJSON[graphs_,numAssociation_Association,opts:OptionsPattern[]]:=Block[{mygraphs=graphs,inMom,outMom,extAsso,cutAsso,llAsso,lnAsso,nameAsso,diagCount=0,procName,exportDir,pLong,fullAsso,diagramwise,processAsso,expAsso,evaluateNumerator,numeratorAsso},
evaluateNumerator[numerator_,numRepl_]:=Block[{numericNum,SP,vector},
  SP[x_List,y_List]:=x[[1]]*y[[1]]-x[[2;;]].y[[2;;]];
  (* vectors are assumed to be covariant *)
  vector[x_List,y_Integer]:=If[y==0,x[[1]],-x[[y+1]] ];
  numericNum=(numerator//. numRepl );
  If[Length@(Variables@numericNum)!=0,
  Print["Error: The numerator coefficients: "<>ToString[Variables@numericNum,InputForm]<>" have no numeric value!"];
  Abort[];  
  ];
  (* short vs long export format *)
  If[!ListQ[numericNum[[1]]],
    numericNum=((numericNum)/.x_/;NumericQ[x]:>ImportString[ExportString[(ReIm@x),"Real64"],"Real64"]),
    numericNum[[All,1]]=(numericNum[[All,1]]/.x_/;NumericQ[x]:>ImportString[ExportString[(ReIm@x),"Real64"],"Real64"]);    
  ];
  numericNum
];
diagramwise=OptionValue[exportDiagramwise];
procName=OptionValue[processName];
exportDir=OptionValue[exportDirectory];
If[ContainsAny[KeyExistsQ[#,"loopLines"]&/@(Flatten@{graphs}),{False}]||ContainsAny[KeyExistsQ[#,"cutStructure"]&/@(Flatten@{graphs}),{False}],
	mygraphs=Flatten@{getCutStructure[graphs]},
	mygraphs=Flatten@{graphs}
];
(* loop over graphs *)
processAsso=Table[
	diagCount+=1;
	(* check input *)
	inMom=Cases[mygraph[["momentumMap"]],_in,Infinity]//Union;
	outMom=Cases[mygraph[["momentumMap"]],_out,Infinity]//Union;
	(* check if stuff is numeric *)
	If[MemberQ[NumericQ/@(Flatten@(inMom /. in[x_]:>x /. numAssociation)),False],
		Print["Error: association is missing numeric value for incoming momenta"];
		Abort[];
	];
	If[MemberQ[NumericQ/@(Flatten@(outMom/. out[x_]:>-x /. numAssociation)),False],
		Print["Error: association is missing numeric value for outgoing momenta"];
		Abort[];
	];
(* check momentum conservation *)
	If[Abs@(Total@((Total@inMom-Total@outMom) /. f_[x_]:>x /. numAssociation))/Abs@(Total@((Total@inMom+Total@outMom+10^-16) /. f_[x_]:>x /. numAssociation))>10^-10,
		Print["Error: rel. error of momentum conservation is worse than 10^(-10)"];
		Abort[]
	];
(* check mass replacements *)
	If[MemberQ[NumericQ/@(Flatten@(mygraph[["loopLines",All,"propagators",All,"mass"]]/.numAssociation)),False],
		Print["No numerical assignment for :"<>ToString[(Flatten@(mygraph[["loopLines",All,"propagators",All,"mass"]]/.numAssociation))/;x_/;NumericQ[x]:>Nothing//Union]];
	];
(* replacement for yaml *)
	extAsso=<|"external_kinematics"->Join[inMom,outMom]/. in[x_]:>x /. out[x_]:>-x/.numAssociation|>;
	cutAsso=<|"ltd_cut_structure"->mygraph[["cutStructure"]]|>;
	(* all the BS is because mathematica exports 0.e-31, which is not readable by yaml *)
	llAsso=
		Map[Evaluate,<|"loop_lines"->Table[
			<|<|"end_node"->69|>,
			<|"propagators"->KeyMap[Replace["mass"->"m_squared"]]/@(l[["propagators"]]/.mass[x_]:>mass[x]^2/.numAssociation /. {x_,y__}/;NumericQ[x]:>ImportString[ExportString[{x,y},"Real64"],"Real64"])|>,
			l[[{"signature"}]],
			<|"start_node"->70|>
			|>,{l,(mygraph[["loopLines"]])}]|>,5];
	lnAsso=Map[Evaluate,<|"n_loops"->Length@mygraph[["loopMomenta"]]|>];
	If[KeyExistsQ[mygraph,"name"],
		nameAsso=<|"name"->mygraph[["name"]]|>,
		nameAsso=<|"name"->procName<>"diag_"<>ToString[diagCount]|>
	];
	If[KeyExistsQ[mygraph,"maximum_ratio_expansion_threshold"],
		expAsso=<|"maximum_ratio_expansion_threshold"->mygraph[["maximum_ratio_expansion_threshold"]]|>,
		expAsso=<|"maximum_ratio_expansion_threshold"->-1|>
	];
	fullAsso=extAsso;
	If[OptionValue[writeNumerator],
	numeratorAsso=<|"numerator_tensor_coefficients"->evaluateNumerator[mygraph[["symmetrizedExpandedTensorCoeff"]],numAssociation]|>;
	numeratorAsso=Evaluate/@numeratorAsso,
	numeratorAsso=<|"numerator_tensor_coefficients"->{{0,0}}|>
	];
	
	Do[fullAsso=Append[fullAsso,asso],{asso,{llAsso,cutAsso,expAsso,lnAsso,numeratorAsso,nameAsso}}];
	If[diagramwise==True,
	If[DirectoryQ[exportDir],
		Export[exportDir<>nameAsso[["name"]]<>".json",fullAsso,"JSON"],
		Print["Couldn't find exportDirectory. Export to standard location."];
		Export["./"<>nameAsso[["name"]]<>".json",fullAsso,"JSON"]
	];
	];
	fullAsso
	,{mygraph,mygraphs}];
	If[DirectoryQ[exportDir],
		Export[exportDir<>"allDiags"<>procName<>".json",processAsso,"JSON"],
		Print["Couldn't find exportDirectory. Export to standard location."];
		Export["./allDiags"<>procName<>".json",processAsso,"JSON"]
	];
	processAsso
]


PrintTemporary["load gamma traces"]

$fileLocation=FileNameJoin[{DirectoryName @ $InputFileName, ""}];

FromFormToMathe[filename_]:=Module[{imp},
		imp=StringJoin[ReadList[filename,String]];
		imp=StringReplace[imp,{"\\"->"","\n"->"","\r" -> "", "("->"[", ")"->"]"," "->"" ,"d_"->"g","mu"->"traceInd"}];
		ToExpression[StringSplit[imp,";"]]
	];





spinorSum[expr_]:=Block[{mom,mass,delta,dummyLong,dummyLong2,lVec,sInd1,sInd2,res},
    gamma[spinorUbar[{mom_,mass_}],dummyLong___,spinorU[{mom_,mass_}]]=gammaTrace[lVec[mom],dummyLong]+mass gammaTrace[dummyLong];
    gamma[spinorVbar[{mom_,mass_}],dummyLong___,spinorV[{mom_,mass_}]]=gammaTrace[lVec[mom],dummyLong]-mass gammaTrace[dummyLong];
    
    gamma[dummyLong__,spinorU[{mom_,mass_}]]gamma[spinorUbar[{mom_,mass_}],dummyLong2__]^=gamma[dummyLong,lVec[mom],dummyLong2]+mass gamma[dummyLong,dummyLong2];
    gamma[dummyLong__,spinorV[{mom_,mass_}]]gamma[spinorVbar[{mom_,mass_}],dummyLong2__]^=gamma[dummyLong,lVec[mom],dummyLong2]-mass gamma[dummyLong,dummyLong2];
    spinorV[{mom_,mass_},sInd1_]spinorVbar[{mom_,mass_},sInd2_]^=gamma[sInd1,lVec[mom],sInd2]-mass delta[sInd1,sInd2] ;
	spinorU[{mom_,mass_},sInd1_]spinorUbar[{mom_,mass_},sInd2_]^=gamma[sInd1,lVec[mom],sInd2]+ mass delta[sInd1,sInd2];
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
gamma[spinorUbar[{mom_,mass_}],dummyLong1___,dummyShort_,lVec[mom_],dummyLong2__]:=-gamma[spinorUbar[{mom,mass}],dummyLong1,lVec[mom],dummyShort,dummyLong2]+antiCom[dummyShort,lVec[mom]] gamma[spinorUbar[{mom,mass}],dummyLong1,dummyLong2];
gamma[spinorUbar[{mom_,mass_}],lVec[mom_],dummyLong1__]:=mass gamma[spinorUbar[{mom,mass}],dummyLong1];
gamma[spinorVbar[{mom_,mass_}],dummyLong1___,dummyShort_,lVec[mom_],dummyLong2__]:=-gamma[spinorVbar[{mom,mass}],dummyLong1,lVec[mom],dummyShort,dummyLong2]+antiCom[dummyShort,lVec[mom]] gamma[spinorVbar[{mom,mass}],dummyLong1,dummyLong2];
gamma[spinorVbar[{mom_,mass_}],mom_,dummyLong1__]:=-mass gamma[spinorVbar[{mom,mass}],dummyLong1];
res1=expr1//Expand
    ];
permuteRight[expr1_]:=Block[{gamma,spinorU,spinorV,res1,lVec},
gamma[dummyLong1__,lVec[mom_],dummyShort_,dummyLong2___,spinorU[{mom_,mass_}]]:=-gamma[dummyLong1,dummyShort,lVec[mom],dummyLong2,spinorU[{mom,mass}]]+antiCom[dummyShort,lVec[mom]] gamma[dummyLong1,dummyLong2,spinorU[{mom,mass}]];
gamma[dummyLong1__,lVec[mom_],spinorU[{mom_,mass_}]]:=mass gamma[dummyLong1,spinorU[{mom,mass}]];
gamma[dummyLong1__,lVec[mom_],dummyShort_,dummyLong2___,spinorV[{mom_,mass_}]]:=-gamma[dummyLong1,dummyShort,lVec[mom],dummyLong2,spinorV[{mom,mass}]]+antiCom[dummyShort,lVec[mom]] gamma[dummyLong1,dummyLong2,spinorV[{mom,mass}]];
gamma[dummyLong1__,lVec[mom_],spinorV[{mom_,mass_}]]:=-mass gamma[dummyLong1,spinorV[{mom,mass}]];
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


Protect[additionalRules,traceFileLocation,spinChainSimplify,polarizationSum];
Options[contractLorentz] = {additionalRules -> {1->1}, traceFileLocation->$fileLocation,spinChainSimplify->False,polarizationSum->False }


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
resultContracted //. addRules
];



contractColor[expr_]:=Block[{delta,T,resultContracted,f,rulesForContraction,TF,Nc,reorderColor,sunA,sunF},
SetAttributes[delta,Orderless];
reorderColor[myexpr_, pats_List] :=
  Module[{h, rls},
    rls = MapIndexed[x : # :> h[#2, Replace[x, rls, -1]] &, pats];
    HoldForm @@ {myexpr /. rls} //. h[_, x_] :> x
  ];
rulesForContraction=Dispatch@{
	delta[sunA[a_],sunA[b_]]f[c__]/;(!FreeQ[{c},a]&&!StringMatchQ[ToString@a,ToString@b]):>ReplaceAll[f[c],a->b],
	f[sunA[a_],sunA[b_],sunA[c_]]:>(-I/TF Module[{n1,n2,n3},T[sunA[a],sunF[n1],sunF[n2]]T[sunA[b],sunF[n2],sunF[n3]]T[sunA[c],sunF[n3],sunF[n1]]-T[sunA[a],sunF[n1],sunF[n2]]T[sunA[c],sunF[n2],sunF[n3]]T[sunA[b],sunF[n3],sunF[n1]]]), (* TF=1/2 *)
	delta[sunA[a_],sunA[b_]]T[c__]/;(!FreeQ[{c},a]&&!StringMatchQ[ToString@a,ToString@b]):>ReplaceAll[T[c],a->b],
	Times[T[a_,i_,j_],T[a_,k_,l_]]:>TF(delta[i,l]delta[j,k]-Power[Nc,-1]delta[i,j]delta[k,l]), (* fierz *)
	T[a_,i_,i_]:>0, (* traceless *)
	T[a_,l_,m_]T[b_,m_,l_]:>TF delta[a,b], (* single trace *)
	delta[a_,b_]delta[b_,c_]/;!StringMatchQ[ToString@a,ToString@c]:>delta[a,c],
	Power[delta[sunA[x_],sunA[y_]],2]/;!StringMatchQ[ToString@x,ToString@y]:>Power[Nc,2]-1,(* adjoint representation *)
	delta[sunA[x_],sunA[x_]]:>Power[Nc,2]-1, (* adjoint *)
	Power[delta[sunF[x_],sunF[y_]],2]/;!StringMatchQ[ToString@x,ToString@y]:>Nc, (* fundamental rep. *)
	delta[sunF[x_],sunF[x_]]:>Nc,(* fundamental *)
	Power[T[b_,c_,d_],2]:>TF(Power[Nc,2]-1),
	Power[f[b_,c_,d_],2]:>2TF(Power[Nc,2]-1)Nc

};

resultContracted=FixedPoint[ ReleaseHold@(ReplaceRepeated[reorderColor[Expand[#],{delta,T,f}],rulesForContraction])&,expr];
resultContracted //.Times[T[a_,i_,j_],T[b_,j_,l_]]/;(!StringMatchQ[ToString@a,ToString@b]&&!StringMatchQ[ToString@i,ToString@l]):>T[a,b,i,l]
];



Protect[consistencyCheckLevel]
Options[extractTensCoeff]={consistencyCheckLevel->1}
extractTensCoeff[graphs_List,opts:OptionsPattern[]]:=extractTensCoeff[#,opts]&/@graphs;
extractTensCoeff[graph_,opts:OptionsPattern[]]:=Block[
	{amp,loopMom,extVect,maxRank,epsK1K1,SP,mySP,tensDecomp,indexSpace,allStruc,allStrucReplRule,allStrucReplRuleVals,res,dummyIndex=Unique[alpha],cleanUpRules
	,tensIndexSet,tensIndices,tensFinal,tensTMP
	,resTMP,replRuleTensExtraction,tens,reorderTensors,repl1,repl2,repl3,tensFinalCheck,finalRes,lMomInd,rewriteRule,ampMod,consistencyCheck=OptionValue[consistencyCheckLevel],$randomDummy},
SetAttributes[SP,Orderless];
loopMom=graph[["loopMomenta"]];
amp=$randomDummy*graph[["numerator"]];
ampMod=amp//. Thread[Rule[lVec/@loopMom,loopMom]];

extVect=Complement[(graph[["momentumMap"]]/. out[x_]:>x /. in[x_]:>x),loopMom];
tensIndexSet=Table[ToExpression@("lInd"<>ToString@lCount),{lCount,Length@loopMom}];


If[amp==0,
	Return[Append[graph,"analyticTensorCoeff"->{{ConstantArray[0,Length@loopMom],0}}]]
];
(* determine max tensorrank *)
maxRank=Max@@(Exponent[(ampMod //. f_[x__]/;!FreeQ[{x},Alternatives@@loopMom]&&StringMatchQ[ToString@f,Alternatives["vector","gamma","SP"]]:>epsK1K1^Count[{x},Alternatives@@loopMom] Unique[f]),epsK1K1]);

tensDecomp=(ampMod//Expand)/. SP[a_,b_]^pow_/;(MemberQ[loopMom,a]||MemberQ[loopMom,b]):>mySP[{a,b},pow]; (* scalarproducts can be to higher powers *)
(* find all structures *)
allStruc=Join[
	Cases[tensDecomp,f_[x__]/;!FreeQ[{x},Alternatives@@loopMom]&&StringMatchQ[ToString@f,Alternatives["vector","gamma","SP"]]:>{f[x],Count[{x} ,Alternatives@@loopMom]},Infinity]//Union,
	Cases[tensDecomp, mySP[x_,pow_]:>{mySP[x,pow],Count[x,Alternatives@@loopMom]*pow} ,Infinity]//Union
	];
(* generate consecutive indexing *)
(* allStruc[[1,-1]] goes to  [startindex(i),endIndex(i)] where startindex(i+1)=endindex(i)+1 *)
allStruc[[1,-1]]={2,allStruc[[1,-1]]+1};
Do[
 allStruc[[countStruc,-1]]={allStruc[[countStruc-1,-1,-1]]+1,allStruc[[countStruc-1,-1,-1]]+allStruc[[countStruc,-1]]}
	,{countStruc,2,Length@allStruc}];
	
	
(* define unique indices *)
allStrucReplRuleVals=Table[
	res={struc[[1]],{}};
	Do[
		res={ReplacePart[#,FirstPosition[#,x_/;MemberQ[loopMom,x]]:>dummyIndex[index]]&@(res[[1]] ),Join[res[[2]],{vector[Extract[#,FirstPosition[#,x_/;MemberQ[loopMom,x]]]&@(res[[1]] ),dummyIndex[index]]}]}
	,{index,Range[struc[[2,1]],struc[[2,-1]]]}
	];
	Times@@(Flatten@res)
	,{struc,(allStruc/. mySP[x_,pow_]:>ConstantArray[SP@@x,pow])}];
	
	
cleanUpRules=Dispatch@{SP[x_,y_]/;(!MemberQ[Join[extVect,loopMom],x]&&!MemberQ[Join[extVect,loopMom],y]):>g[x,y],
						 SP[x_,y_]/;(!MemberQ[Join[extVect,loopMom],x]&&MemberQ[Join[extVect,loopMom],y]):>vector[y,x],
						 SP[x_,y_]/;(MemberQ[Join[extVect,loopMom],x]&&!MemberQ[Join[extVect,loopMom],y]):>vector[x,y]
						 };	
allStrucReplRuleVals=allStrucReplRuleVals //.cleanUpRules ;
allStrucReplRule=Thread@(Rule[allStruc[[All,1]],allStrucReplRuleVals]);
tensDecomp=(tensDecomp//. allStrucReplRule);

(* consistency check --> new expression matches input: Extraction of loop-momenta worked *)
If[consistencyCheck>1,
	Print["check I started"];
	If[Simplify@((amp-contractLorentz[tensDecomp]) /. Thread[Rule[Join[loopMom,extVect],Array[Prime,Length@(Join[loopMom,extVect])]]])=!=0
		,
		Print["Error in tensor extraction: Check Input"]; Abort[]
		,
		Print["check I passed"]
	 ]
];

tensDecomp=Table[SeriesCoefficient[(tensDecomp /. vector[x_,ind_]/;MemberQ[loopMom,x]:>epsK1K1 vector[x,ind]),{epsK1K1,0,rank}],{rank,0,maxRank}];
(* consistency check --> new expression matches input: rank extraction worked *)
If[consistencyCheck>1,
	Print["check II started"];
	If[Simplify@((amp-Total@(contractLorentz[tensDecomp]) /. Thread[Rule[Join[loopMom,extVect],Array[Prime,Length@(Join[loopMom,extVect])]]]))=!=0
		,
		Print["Error in tensor extraction: Check Input"]; Abort[]
		,
		Print["check II passed"]
	 ]
];
(* redefine mathematica default ordering in times *)
reorderTensors[myexpr_, pats_List] :=
  Module[{h, rls},
    rls = MapIndexed[x : # :> h[#2, Replace[x, rls, -1]] &, pats];
    HoldForm @@ {myexpr /. rls} //. h[_, x_] :> x
  ];

tensFinal=Table[
	resTMP=tensDecomp[[rankCount+1]];
	Do[
	tensIndices=Table[tensIndexSet[[loopCount]][ind],{ind,1,rankCount}];
	(*works since we defined unique indices*)
	tensTMP=Cases[tensDecomp[[rankCount+1]],vector[loopMom[[loopCount]],ind_]:>ind,Infinity]//Union;
	(*Print[tensTMP];*)
		replRuleTensExtraction=Dispatch@{
		tens[x_,y_]tens[z_,v_]:>tens[Join[x,z],Join[y,v]]		 	
	};
	repl1=f_[x___]/;(MemberQ[{vector,g,gamma},Head@(f[x])]&&!FreeQ[Flatten@{x},Alternatives@@tensTMP]):>tens[{f[x]}, Cases[{x},Alternatives@@tensTMP,Infinity]];
	repl2=tens[{start___,tens[a_,x_],rest___} ,y_]:> tens[Flatten@Join[{start},{a},{rest}], Join[x,y]];
	repl3 = tens[x_,y_]/;!FreeQ[y,Alternatives@@tensTMP] :> ReplaceAll[tens[x, y],Thread@(Rule[Cases[y,Alternatives@@tensTMP,Infinity],Flatten@{tensIndices[[(*rankCount+1,*)1;;Length@Cases[y, Alternatives@@tensTMP,Infinity]]]} ] )];
	resTMP=(ReleaseHold@(ReplaceRepeated[reorderTensors[Expand[(# //.vector[loopMom[[loopCount]],ind_]-> 1 /. repl1//. repl2)],{tens}],replRuleTensExtraction])&@((resTMP)))/.repl3;  
	,{loopCount,Length@loopMom}];
	resTMP
,{rankCount,0,maxRank}];

(* order entries by loop-momenta tensors *)
tensFinal=Flatten[Table[
	Table[
		{tensCoeff,(tensFinal[[rankCount+1]] //.tens[x_,ind_]/;(Count[ind,#]&/@(Blank/@(tensIndexSet[[;;Length@loopMom]]))!=tensCoeff):>0 ) }
		,{tensCoeff,Sort[(Flatten[#,1]&@(Permutations/@(PadRight[#,Length@loopMom]&/@IntegerPartitions[rankCount,Length@loopMom])))]//Reverse}]	
,{rankCount,0,maxRank}],1];
tensFinal=tensFinal //. tens[x_,y_]:>Times@@x;

tensFinal=tensFinal //. {loopList_List,{x_}}/;(Length@loopList==Length@loopMom && FreeQ[loopList,Alternatives@@extVect]):>{loopList,x};
 tensFinal=DeleteCases[tensFinal,{x_List,0},1];
 (* consistency check --> final expression matches input *)
(* final check*)
If[consistencyCheck>=1,
	Print["Final check started"];
	(* replace integers by loop-momenta *)
	tensFinalCheck=(tensFinal//. {indX_List,x_}/;AllTrue[indX,IntegerQ]:>({loopMom^2 loopMom^indX,x}//.mom_^pow_/;MemberQ[loopMom,mom]:>(Times@@(Array[Evaluate@vector[mom,Extract[tensIndexSet,Position[loopMom, mom]][[1]][#-2]]&,pow] ))))//.vector[ll_,_[num_]]/;num <= 0:>1;
	(* restore input structure *)
	tensFinalCheck=contractLorentz[(tensFinalCheck //. {loopList_List,x_}/;(Length@loopList==Length@loopMom && FreeQ[loopList,Alternatives@@extVect]):>(Times@@(loopList)*x))];
	tensFinalCheck=(Total@tensFinalCheck);
	
	If[Simplify@((amp-tensFinalCheck) /. Thread[Rule[Join[loopMom,extVect],Array[Prime,Length@(Join[loopMom,extVect])]]])=!=0
	,
	Print["Error in final check: likely a bug :("]; Abort[],
	Print["Final check passed"]
	]
	];
	
	(* replace lIndX[Y] by lMomInd[x][y]... x: loop-momentum number, y: lorentz-index: mu[y]...  ugly but was only fix which came to mind so \.af\_(\:30c4)_/\.af  *)
	rewriteRule=Thread[RuleDelayed[Compose[#[[1]],#[[2]]]&/@(Transpose@{tensIndexSet,ConstantArray[Pattern[x,Blank[]],Length@tensIndexSet]}),Evaluate[Compose[#[[1]],#[[2]]]&/@(Transpose@{Array[lMomInd,Length@tensIndexSet],ConstantArray[x,Length@tensIndexSet]})]]];
	
	tensFinal=tensFinal/.rewriteRule /. $randomDummy->1;
	finalRes=Append[graph,"analyticTensorCoeff"->tensFinal]
]



Protect[outputFormat];
Options[getSymCoeff] ={ outputFormat->"long"};

getSymCoeff[graphs_,opts:OptionsPattern[]]:=Block[
{tensDecompNum,resTmp,res,indSet,g,rank,replV,replG,format,indShift,replGamma,fullResult,resTmpTmp},
    format=OptionValue[outputFormat];	
	g[int1_Integer,int2_Integer]:=If[int1!=int2,0, If[int1==0,1,-1]];
	fullResult=Table[
    If[KeyExistsQ[graph,"analyticTensorCoeff"],
      tensDecompNum=graph[["analyticTensorCoeff"]]
      ,
      tensDecompNum=extractTensCoeff[graph];
      tensDecompNum=tensDecompNum[["analyticTensorCoeff"]]
    ];
    If[format==="long",
      rank=tensDecompNum[[All,1]];
      (* generate complete rank upto max rank *)
      rank=Flatten[Table[Sort@Flatten[Permutations/@IntegerPartitions[cR,{Length@rank[[-1]]},Range[0,cR]],1],{cR,0,Total@rank[[-1]]}],1]
      ,
      (* only non-vanishing coefficients*)
      rank=tensDecompNum[[All,1]]
    ];
	
   res=
   Flatten[
   Table[
    (* translate given rank to all possible (symmetrized sets) of lorentz momenta *)
	indSet=Distribute[(rank[[rankCount]]/. x_/;NumericQ[x]:>(Permutations/@(Sort@Tuples[Range[0,3],{x}]))),List];
	
	If[MemberQ[tensDecompNum[[All,1]],rank[[rankCount]]],
		resTmp=Extract[tensDecompNum[[All,2]],Position[tensDecompNum[[All,1]],rank[[rankCount]]]];
		resTmp=resTmp[[1]]
		,
		resTmp=0	
	];
	indShift=Table[ConstantArray[(llCount-1)*4,(rank[[rankCount,llCount]])],{llCount,Length@rank[[rankCount]]}];
	(* loop over all symmetrized sets *)
	Table[
	    replV=vector[x_,lMomInd[loopMomNum_][loopIndex_]]:>Sum[vector[x,sInd],{sInd,myInd[[loopMomNum,All,loopIndex]]}];
		replG=g[lMomInd[loopMomNum_][loopIndex_],lMomInd[loopMomNum2_][loopIndex2_]]:>Sum[g[s1,s2],{s1,myInd[[loopMomNum2,All,loopIndex2]]},{s2,myInd[[loopMomNum2,All,loopIndex2]]}];
		replGamma=gamma[x___,lMomInd[loopMomNum_][loopIndex_],y___]:>Sum[gamma[x,sInd,y],{sInd,myInd[[loopMomNum,All,loopIndex]]}];
	    resTmpTmp=resTmp //. replV //. replG //. replGamma;
	    {resTmpTmp,Flatten@(myInd[[All,1]]+indShift)}
	,{myInd,indSet}]
  ,{rankCount,Length@rank}]
  ,1];

(*res=res //. gamma[x__,y_,z__]/;(!Head[y]===lVec&&!NumericQ[y])\[RuleDelayed]Sum[(gamma[x,y,z]/. y\[Rule]lInd),{lInd,0,3}];*)

    If[format==="long",
      res=res[[All,1]]
      ,
      (* only non-vanishing coefficients*)
      res=DeleteCases[res,{0,x_List},1]
    ];
    
Append[graph,"symmetrizedExpandedTensorCoeff"->res]
,{graph,Flatten@{graphs}}];

If[Length@fullResult==1,
fullResult[[1]],
fullResult
]
]



Options[processNumerator]={algebraOptions->{additionalRules -> {1->1}, traceFileLocation->$fileLocation,spinChainSimplify->False,polarizationSum->False },extractTensors->False, symCoefficients->True,coeffFormat->"long"};
processNumerator[graphs_,pathToFeynmanRules_,opts:OptionsPattern[]]:=Block[{myGraphs=Flatten@{graphs},loopMom},
If[FileExistsQ[pathToFeynmanRules],
Get[pathToFeynmanRules],
Print["Error: Could not find Feynman rules in: \""<>pathToFeynmanRules<>"\""];
Abort[]
];
loopMom=myGraphs[[1,"loopMomenta"]];
myGraphs=Map[Evaluate,#]&/@myGraphs//. scalarProp[x_,y_]/;!FreeQ[x,Alternatives@@loopMom]:>1 //. scalarProp[x_,y_]/;FreeQ[x,Alternatives@@loopMom]:>1 /(SP[x,x]-y^2);
myGraphs=Map[Evaluate,#]&/@myGraphs;
myGraphs[[All,"numerator"]]=contractColor[myGraphs[[All,"numerator"]]];
myGraphs[[All,"numerator"]]=contractLorentz[myGraphs[[All,"numerator"]],OptionValue[algebraOptions]];
myGraphs=DeleteCases[myGraphs,x_/;x[["numerator"]]==0,1];
If[OptionValue[extractTensors]==True,
myGraphs=extractTensCoeff[myGraphs]
];
If[OptionValue[symCoefficients]==True,
myGraphs=getSymCoeff[myGraphs,outputFormat->OptionValue[coeffFormat]]
];
If[Length@myGraphs==1,
	myGraphs[[1]],
	myGraphs
]
]
