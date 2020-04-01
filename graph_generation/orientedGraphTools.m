(* ::Package:: *)

Print["The documented functions in this package are: \n ?plotGraph \n ?findIsomorphicGraphs \n ?constructCuts
 ----------------------------------------- 
 Needs the package IGraphM which can be downloaded from https://github.com/szhorvat/IGraphM. !!! \n Run: Get[\"https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m\"] for installation"
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


ClearAll[constructCuts]
Protect[NumberFinalStates,DisplayCuts,DetectSelfEnergy];
Options[constructCuts] ={ NumberFinalStates->All,DisplayCuts->True,DetectSelfEnergy-> True};
constructCuts[graph_,opts:OptionsPattern[]]:=Block[{finalStateNum,detectSE,displayCuts,particles,momenta, trees,edges,cutEdges,forests,external,vertices,tmpDiag,cutTag,loopNum,cutGraphs,seInfo,props,cutProps,cutPos,prop,seEdges,probEdge,selfEnergy,subGEdges,sePostion,sePos,findAllSpanningTrees},
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
finalStateNum=OptionValue[NumberFinalStates];
displayCuts=OptionValue[DisplayCuts];
detectSE=OptionValue[DetectSelfEnergy];

edges=graph[["edges"]]/. Rule[a_,b_]:>UndirectedEdge[a,b] /. DirectedEdge[a_,b_]:>UndirectedEdge[a,b]/. TwoWayRule[a_,b_]:>UndirectedEdge[a,b];
momenta=graph[["momentumMap"]];
particles=graph[["particleType"]];
external=Cases[edges,UndirectedEdge[a_,b_]/;!FreeQ[a,out]||!FreeQ[a,in]||!FreeQ[b,out]||!FreeQ[b,in],Infinity]//Union;
vertices=VertexList@edges;
loopNum=Length@(Complement[Variables@momenta,Variables@({Extract[momenta,Position[edges,Alternatives@@external]],Extract[momenta,Position[edges,Alternatives@@external]]/. in[x_]:>x /. out[x_]:>x} )]);


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
			{shifts[[countShifts,2,isoCount]],FindPermutation[graphsEdgeColored[[shifts[[countShifts,1,1]]]],graphsEdgeColored[[shifts[[countShifts,1,2]]]]/.shifts[[countShifts,2,isoCount]]]}
		,{isoCount,Length@shifts[[countShifts,2]]}]}]
	,{countShifts,Length@shifts}];
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
	shiftsAndPermutations
];


(* GraphComputation`GraphPlotLegacy is a hack, since mathematica updated GraphPlot *)
Protect[edgeLabels,plotSize]
Options[plotGraph] ={ edgeLabels->{},plotSize->Scaled[0.5]};
plotGraph[graph_,opt:OptionsPattern[]]:=Block[{mygraph,imSize=OptionValue[plotSize]},
		
		mygraph=Transpose@(Values@(graph[[Join[{"edges"},Flatten@{OptionValue[edgeLabels]}]]]))/. TwoWayRule->Rule /. UndirectedEdge->Rule /. DirectedEdge->Rule /. {Rule[a_,b_],c__}:>{Rule[a,b],Flatten@{c}};
		
	Print[GraphComputation`GraphPlotLegacy[mygraph,Method->{"SpringElectricalEmbedding", "RepulsiveForcePower" -> -1.55},EdgeRenderingFunction->({{RGBColor[0.22,0.34,0.63],Arrowheads[0.015],Arrow[#]},If[#3=!=None,Text[ToString[#3],Mean[#1],Background->White],{}]}&),MultiedgeStyle->.3,SelfLoopStyle->All,ImageSize->imSize]]
];
