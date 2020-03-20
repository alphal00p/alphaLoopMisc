(* ::Package:: *)


Print["The documented functions in this package are: \n ?plotGraph \n ?findIsomorphicGraphs
 ----------------------------------------- 
 Needs the package IGraphM which can be downloaded from https://github.com/szhorvat/IGraphM. !!! \n Run: Get[\"https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m\"] for installation"
];
Needs["IGraphM`"];
findIsomorphicGraphs::usage="Finds isomorphisms between graphs.
	Input: List of graphs
	Represent graph by: {{edge1: (Rule[a,b] for directed edge, TwoWayRule[a,b] for undirected edge), tag1:List},...,{edgeN, tagN:List}} 
	No tagging: empty list,
	Output: List where every entry has the form
	{
		{Position of graph to be mapped, Positon of graph it maps to},
		{
		{isomorphism, permutation of edges, orientation changes: +1 = no orientation change, -1: orientation change}
		}
	}
"
plotGraph::usage "plots graph
	Input plotGraph[graph, image size, default: Scaled[0.5]]
	Represent graph by: {{edge1: (Rule[a,b] for directed edge, TwoWayRule[a,b] for undirected edge),edgelabel},...,{edgeN,edgelabelN}} "



(* represent graph by: {{edge1: (Rule[a,b] for directed edge, TwoWayRule[a,b] for undirected edge),  ,tag1:List},...,{edgeN,tagN:List}} *)

findIsomorphicGraphs[graphs_:List]:=Block[{tagging,taggingToNum,$uniqueTag, graphsEdgeColored,graphsEdgeColoredSimple,isIsomorph,twr,$myAssoc,isomorphisms,shifts,g1,g2,shiftsAndPermutations,orientation,mappedGraph,graphID1,graphID2,permutation,mapping},
	SetAttributes[twr,Orderless];
	(* taggings become colors on edge-colored graphs, multiple taggings will be translated to single prime which becomes the color of the edge *)
	tagging=Variables@(graphs[[All,All,2]]);
	If[Length@tagging!=0,
		tagging=tagging,
		tagging=$uniqueTag
	];
	(* Quick Input check of input *)
	If[MemberQ[ListQ/@(Flatten[graphs[[All,All,2]],1]),False],
	Print["Error: tags of edges have to be given as a list. No tagging: {}"];
	Abort[]		
	];
	
(* start prime range late, so that translation of multi-edge graphs is still unique *)
	taggingToNum=Thread@(Rule[tagging,Prime[Range[100,Length@tagging+99]]]);
(*graphIsomorph=Table[IGIsomorphicQ[Graph@graph[i],Graph@graph[j]],{i,1,Length[momenta]},{j,1,Length[momenta]}]//LowerTriangularize;
isomorphisms=(Position[graphIsomorph,True]/. {i_,i_}-> Nothing);
*)
(* translate to undirected graphs and apply tagging *)
	graphsEdgeColored=graphs /. taggingToNum /. {Rule[a_,b_],x_List}:>{twr[a,b],Times@@x} /. {TwoWayRule[a_,b_],x_List}:>{twr[a,b],Times@@x};

(* get isomorphism for graphs which are isomorphic; IGVF2 needed for multi-edge graphs: 
	There is no algorithm working for multi-edge graphs directly. Therefore translate into edge-colored simple graphs! The color is the multiplicity of the edge e.g.: g1Multiedge={{1\[TwoWayRule]4,5},1\[TwoWayRule]4,5}} --> g1ColouredSimple={1\[TwoWayRule]4\[Rule]2*5}*)
	graphsEdgeColoredSimple=(Counts/@(Sort/@graphsEdgeColored)) /.Association:>$myAssoc/.Rule[List[twr[a_,b_],color1_],color2_]:> Rule[twr[a,b],color1*color2] /. twr->TwoWayRule /. $myAssoc:>Association;
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
			mappedGraph=Permute[graphs[[graphID1]],permutation]/.mapping /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b});
			orientation=(Keys[mappedGraph[[All,1]]]-Values[mappedGraph[[All,1]]])/(Keys[(graphs[[graphID2,All,1]] /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b}))]-Values[(graphs[[graphID2,All,1]] /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b}))])//Simplify;
			shiftsAndPermutations[[shiftCount,2,isoCount]]=Join[shiftsAndPermutations[[shiftCount,2,isoCount]],{orientation}];
		,{isoCount,Length@shiftsAndPermutations[[shiftCount,2]]}];
	,{shiftCount,Length@shiftsAndPermutations}];
	shiftsAndPermutations
];
(* GraphComputation`GraphPlotLegacy is a hack, since mathematica updated GraphPlot *)
plotGraph[graph_,imSize_:Scaled[0.5]]:=Block[{mygraph=graph},
	If[Union@(Length/@mygraph)!={2},
		Print["No or incomplete labeling. Automatic labeling activated"];
		mygraph=Transpose@{mygraph,Range[Length@mygraph]}	
	];
	Print[GraphComputation`GraphPlotLegacy[mygraph,Method->{"SpringElectricalEmbedding", "RepulsiveForcePower" -> -1.55},EdgeRenderingFunction->({{RGBColor[0.22,0.34,0.63],Arrowheads[0.015],Arrow[#]},If[#3=!=None,Text[#3,Mean[#1],Background->White],{}]}&),MultiedgeStyle->.3,SelfLoopStyle->All,ImageSize->imSize]]
];
