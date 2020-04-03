(* ::Package:: *)

Print["The documented functions in this package are: \n ?plotGraph \n ?findIsomorphicGraphs \n ?constructCuts \n ?importGraphs \n ?getLoopLines \n ?getCutStructure \n ?writeMinimalJSON
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
importGraphs::usage="Import QGraf graphs generated with orientedGraphs.sty
			Input: 
				\"PathToQGrafOutput\"
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
			Options: {processName\[Rule] \"NameOfProcess\"(Default: \"\"), exportDirectory\[Rule] \"PathToExportDir/\" (Default: \"./\")}
		Output:
			None
"


ClearAll[getLoopLines]
getLoopLines[graphs_]:=Block[{eMom,lMom,lLines,lProp,orientation,res,mass},
res=Table[
If[!KeyExistsQ[graph,"loopMomenta"],
	Print["Error: Graphs have no entry \"loopMomenta\". Import with importGraph[]"];
	Abort[]
];
eMom=Complement[Variables@(Union@Flatten@graph[["momentumMap"]]),graph[["loopMomenta"]]];
lLines=CoefficientArrays[((graph[["momentumMap"]] /. Thread[Rule[eMom,0]]/. 0:>Nothing)^2 //FullSimplify)/.(x_)^2:>x// DeleteDuplicates,graph[["loopMomenta"]]][[2]]//Normal;

lProp=Table[
	orientation=(ReplaceAll[#,x_/;!NumberQ[x]:>0]&/@(FullSimplify@((graph[["momentumMap"]]/. Thread[Rule[eMom,0]])/(Dot[graph[["loopMomenta"]],ll])) ));
	{ll,(graph[["momentumMap"]]orientation)/. 0:>Nothing,mass/@(graph[["particleType"]]Abs[orientation])/. mass[0]:>Nothing} /. Thread[Rule[graph[["loopMomenta"]],0]]
,{ll,lLines}];
Append[graph,<|"loopLines"->lProp|>]
,{graph,Flatten@{graphs}}];
If[Length@res==1,
res=res[[1]]
];
res
]


ClearAll[importGraphs]
importGraphs[file_:String]:=Block[{graphs, loopMom},
graphs=ToExpression["{"<>StringDrop[StringReplace[Import[file,"Text"],{"**NewDiagram**"->",","\n"->""}]<>"}",1]];
loopMom=Complement[Variables@(Union@Flatten@graphs[[All,"momentumMap"]]/.f_[x_]:>x),Union@(Cases[graphs[[All,"momentumMap"]],f_[x_]:>x,Infinity] )]
;
graphs=Append[#,"loopMomenta"->loopMom]&/@graphs
]


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
	shiftsAndPermutations
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


ClearAll[getCutStructure]

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
	ll=gg[["loopLines"]][[All,1]];
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


ClearAll[writeMinimalJSON]
Protect[processName,exportDirectory];
Options[writeMinimalJSON]={processName->"",exportDirectory->"./"}
writeMinimalJSON[graphs_,numAssociation_Association,opts:OptionsPattern[]]:=Block[{mygraphs=graphs,inMom,outMom,extAsso,cutAsso,llAsso,lnAsso,nameAsso,diagCount=0,procName,exportDir,pLong,fullAsso},
procName=OptionValue[processName];
exportDir=OptionValue[exportDirectory];
If[ContainsAny[KeyExistsQ[#,"loopLines"]&/@(Flatten@{graphs}),{False}]||ContainsAny[KeyExistsQ[#,"cutStructure"]&/@(Flatten@{graphs}),{False}],
	mygraphs=Flatten@{getCutStructure[graphs]},
	mygraphs=Flatten@{graphs}
];
(* loop over graphs *)
Table[
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
	If[MemberQ[NumericQ/@(Flatten@(mygraph[["loopLines"]]/.numAssociation)),False],
		Print["No numerical assignement for :"<>ToString[(Flatten@(mygraph[["loopLines"]]/.numAssociation))/;x_/;NumericQ[x]:>Nothing//Union]];
	];
(* replacement for yaml *)
	extAsso=<|"external_kinematics"->Join[inMom,outMom]/. in[x_]:>x /. out[x_]:>-x/.numAssociation|>;
	cutAsso=<|"ltd_cut_structure"->mygraph[["cutStructure"]]|>;
	llAsso=
		<|"looplines"->Table[
			<|<|"end_note"->69|>,
			<|"propagators"->Flatten@Table[
				If[Length@l[[-2,p]]==1,pLong= Flatten@ConstantArray[l[[-2,p]],4],pLong=Flatten@l[[-2,p]]];
				<|<|"m_squared"->l[[-1,p]][[1]]|>,
				<|"q"->pLong|>
				|>
				,{p,Length@l[[2]]}]|>
			,<|"signature"->l[[1]]|>
			,<|"start_note"->70|>
			|>,{l,(mygraph[["loopLines"]]/. mass[x_]:>mass[x]^2/.numAssociation /. x_/;NumericQ[x]:>ImportString[ExportString[x,"Real64"],"Real64"])}]|>;
	lnAsso=Map[Evaluate,<|"n_loops"->Length@mygraph[["loopMomenta"]]|>];
	If[KeyExistsQ[mygraph,"name"],
		nameAsso=<|"name"->mygraph[["name"]]|>,
		nameAsso=<|"name"->procName<>"diag_"<>ToString[diagCount]|>
	];
	fullAsso=extAsso;
	Do[fullAsso=Append[fullAsso,asso],{asso,{llAsso,cutAsso,lnAsso,nameAsso}}];
	If[DirectoryQ[exportDir],
		Export[exportDir<>nameAsso[["name"]]<>".json",fullAsso,"JSON"],
		Print["Couldn't find exportDirectory. Export to standard location."];
		Export["./"<>nameAsso[["name"]]<>".json",fullAsso,"JSON"]
	];
	,{mygraph,mygraphs}];
fullAsso
]
