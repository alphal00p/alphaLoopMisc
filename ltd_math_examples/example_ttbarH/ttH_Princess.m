(* ::Package:: *)

NoteBookDir="/home/hirschva/MG5/alphaLoopMisc/ltd_math_examples/example_ttbarH";
LTDToolsPath="/home/hirschva/MG5/alphaLoopMisc/ltd_math_utils/LTDTools.m";


SetDirectory[NoteBookDir];


Get[LTDToolsPath];


allTTHGraphs=importGraphs["./epem_ttH.qgr",sumIsoGraphs->False];


SelectedDiagram=allTTHGraphs[[65]];


(* ::Subsection::Initialization:: *)
(*(* Building loop momenta basis *)*)


(* ::Input::Initialization:: *)
BuildMomentumBasisFromDefiningPropagators[definingPropagators_]:=Block[
{
iProp,
linearSystem={},
kis={k1,k2,k3,k4},
kisp={k1p,k2p,k3p,k4p}
},
For[iProp=1,iProp<=Length[definingPropagators],iProp++,
AppendTo[linearSystem,Evaluate[kisp[[iProp]]]==SelectedDiagram["momentumMap"][[definingPropagators[[iProp]]]]];
];
Solve[linearSystem,kis][[1]]/.Table[
Evaluate[kisp[[j]]]->Evaluate[kis[[j]]]
,{j,Range[Length[kis]]}]
];


(* ::Input::Initialization:: *)
combinationsOfPropagatorsDefiningTheMomentumBasisOfEachCutkoskyCut=<|
{10,11,12,14,16}->{10,11,12,14}, (* #1 RR x RR *)
{10,11,15} ->{10,11,14,16}, (* #2 VV x B *)
{10,11,16,17}->{10,11,16,14}, (* #3 RV x R *)
{11,12,13,14}->{11,12,13,16}, (* #4 R x RV *)
{11,12,8}->{11,12,14,16}, (* #5 B x VV *)
{11, 13, 15, 16} -> {11,13,15,14}, (* #6 RV x R *)
{11,13,17} ->{11,13,14,16}, (* #7 V x V *)
{11,14,15,16,8}->{11,14,15,16}, (* #8 RR x RR *)
{11,14,17,8}->{11,14,17,16} (* #9 R x RV *)
|>;


(* ::Input::Initialization:: *)
AllSelectedCutkoskyCutsMomentumBases=Association[Table[
key->BuildMomentumBasisFromDefiningPropagators[combinationsOfPropagatorsDefiningTheMomentumBasisOfEachCutkoskyCut[key]],
{key,Keys[combinationsOfPropagatorsDefiningTheMomentumBasisOfEachCutkoskyCut]}]];


(* ::Input::Initialization:: *)
FromRustToQGraf={
k1->k1,
k2->k1+k4-p1-p2,
k3->-k3,
k4->-k2+k4-p1-p2
};


(* ::Input::Initialization:: *)
FromQGrafToRustRotation=(Solve[
{
k1p==k1,
k2p==k1+k4-p1-p2,
k3p==-k3,
k4p==-k2+k4-p1-p2
}
,
{k1,k2,k3,k4}
][[1]])/.{k1p->k1,k2p->k2,k3p->k3,k4p->k4};


(* ::Input::Initialization:: *)
CombineRotations[rot1_,rot2_]:=Block[
{},
{
k1->Simplify[((k1/.rot1)/.rot2)],
k2->Simplify[((k2/.rot1)/.rot2)],
k3->Simplify[((k3/.rot1)/.rot2)],
k4->Simplify[((k4/.rot1)/.rot2)]
}
]


(* ::Input::Initialization:: *)
FinalAllSelectedCutkoskyCutsMomentumBases=Association[Table[
key->CombineRotations[AllSelectedCutkoskyCutsMomentumBases[key],FromQGrafToRustRotation]
,{key,Keys[AllSelectedCutkoskyCutsMomentumBases]}]];


(* ::Text:: *)
(*Consistency check*)


(* ::Input::Initialization:: *)
If[Table[
Table[Simplify[SelectedDiagram["momentumMap"][[iProp]]/.FinalAllSelectedCutkoskyCutsMomentumBases[key]/.FromRustToQGraf],
{iProp,combinationsOfPropagatorsDefiningTheMomentumBasisOfEachCutkoskyCut[key]}
]==={k1,k2,k3,k4}
,{key,Keys[AllSelectedCutkoskyCutsMomentumBases]}
]==Table[True,{D,Range[Length[FinalAllSelectedCutkoskyCutsMomentumBases]]}],Print["Consistency check PASSED!"];,Print["Consistency check FAILED!"];];


(* ::Subsection::Initialization:: *)
(*(*Process numerator coefficients*)*)


Print["Now performing algebra of this diagram..."];


processedSelectedGraph=Import["./processedSelectedGraphNew.m"];
Print["Done!"];


(* ::Subsubsection::Initialization:: *)
(*(*Generate coefficients*)*)


(* ::Input::Initialization:: *)
MyPrefactor = I ge^4 gs^4 yt^2;


(* ::Input::Initialization:: *)
SelectedGraphNumeratorAsList=List@@processedSelectedGraph["numerator"];


(* ::Input::Initialization:: *)
ApplyFunctionToEachTerm[terms_,func_,OptionsPattern[{DEBUG->True,printoutFrequency->1,msgPrefix->"",termIndexOffset->0}]]:=Block[
{i, newTerms={},temp},
For[i=1,i<=Length[terms],i++,
If[And[OptionValue[DEBUG],Mod[i,OptionValue[printoutFrequency]]==0],
NotebookDelete[temp];
temp=PrintTemporary[OptionValue[msgPrefix]<>"Now processing terms >#"<>ToString[i+OptionValue[termIndexOffset]]<>"/"<>ToString[Length[terms]+OptionValue[termIndexOffset]]<>"  "<>ToString[N[i/Length[terms],3]*100.]<>"% "<>"..."]];
AppendTo[newTerms,func[terms[[i]]]];
];
NotebookDelete[temp];
newTerms
]


(* ::Text::Initialization:: *)
(*(*Substitute in numerical values*)*)


(* ::Input::Initialization:: *)
ProcessedSelectedGraphNumeratorAsList=List@@(Total[SelectedGraphNumeratorAsList]/.Dispatch[{mH->125,mT->173,Nc->3, d->4,SP[p1,p2]:>2 10^6}]);


(* ::Input::Initialization:: *)
ProcessedSelectedGraphNumeratorAsList=ApplyFunctionToEachTerm[ProcessedSelectedGraphNumeratorAsList,(Simplify[#/MyPrefactor])&,DEBUG->False,printoutFrequency->1000];


(* ::Input::Initialization:: *)
MyEvaluateNumerator[numerator_,numRepl_]:=Block[{numericNum=numerator,Pair},
  numericNum=numericNum//FCI;
  Pair[Momentum[x_List,d___],Momentum[y_List,d___]]:=x[[1]]*y[[1]]-x[[2;;]].y[[2;;]];
  (* vectors are assumed to be covariant *)
   numericNum=(numerator//. numRepl //.x_List[y_Integer]:>If[y==0,x[[1]],-x[[y+1]] ]);
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


(* ::Input::Initialization:: *)
BuildNumericalSymmetricCoeffs[analyticalExpression_,OptionsPattern[{consistencyCheckLevel->0,ReplacementRules->{}}]]:=Block[
{replRules={
p1->{1000.0,0,0,1000.0},
p2->{1000.0,0,0,-1000.0},
mH->125.0,mT->173.0,Nc->3, d->4
},
coeffs
},
coeffs=MyEvaluateNumerator[
getSymCoeff[<|
"numerator"->(ExpandScalarProduct[analyticalExpression/.OptionValue[ReplacementRules]]),
"momentumMap"->SelectedDiagram["momentumMap"],
"loopMomenta"->SelectedDiagram["loopMomenta"]
|>,outputFormat->"short",consistencyCheckLevel->OptionValue[consistencyCheckLevel]]["symmetrizedExpandedTensorCoeff"],
{
p1->{1000,0,0,1000},
p2->{1000,0,0,-1000},
mH->125,mT->173,Nc->3, d->4
}
];
(*coeffs*)
DeleteCases[coeffs,x_/;And[x[[1]][[1]]==0,x[[1]][[2]]==0]]
]


(* ::Input::Initialization:: *)
SumCoeffs[coeffs_]:=Module[{SummedCoeffs,iCoeff,iCoeffsList, coeffPos},
SummedCoeffs={};
For[iCoeffsList=1,iCoeffsList<=Length[coeffs],iCoeffsList++,
For[iCoeff=1, iCoeff<=Length[coeffs[[iCoeffsList]]],iCoeff++,
coeffPos=FirstPosition[SummedCoeffs,{{_,_},coeffs[[iCoeffsList]][[iCoeff]][[2]]}];
If[coeffPos===Missing["NotFound"],
AppendTo[SummedCoeffs,coeffs[[iCoeffsList]][[iCoeff]]];,
SummedCoeffs[[coeffPos[[1]]]][[1]][[1]]=SummedCoeffs[[coeffPos[[1]]]][[1]][[1]]+coeffs[[iCoeffsList]][[iCoeff]][[1]][[1]];
SummedCoeffs[[coeffPos[[1]]]][[1]][[2]]=SummedCoeffs[[coeffPos[[1]]]][[1]][[2]]+coeffs[[iCoeffsList]][[iCoeff]][[1]][[2]];
];
];
];
SummedCoeffs=DeleteCases[SummedCoeffs,x_/;And[x[[1]][[1]]==0,x[[1]][[2]]==0]];
SortBy[SummedCoeffs,(#[[2]])&]
]


(* ::Text::Initialization:: *)
(*(*Build all coefficients*)*)


(* ::Input::Initialization:: *)
GenerateNumeratorCoefficients[numeratorAsListOfTerms_,OptionsPattern[{
consistencyCheckLevel->0,ReplacementRules->{},DEBUG->True,printoutFrequency->10000,msgPrefix->"",termIndexOffset->0
}]
]:=Block[
{overallCoeffs={},functionToApply},
functionToApply[x_]:=overallCoeffs=SumCoeffs[{overallCoeffs,BuildNumericalSymmetricCoeffs[x,
ReplacementRules->OptionValue[ReplacementRules],consistencyCheckLevel->OptionValue[consistencyCheckLevel]
]}];
ApplyFunctionToEachTerm[numeratorAsListOfTerms,functionToApply,
DEBUG->OptionValue[DEBUG],printoutFrequency->OptionValue[printoutFrequency],msgPrefix->OptionValue[msgPrefix],termIndexOffset->OptionValue[termIndexOffset]
];
overallCoeffs
]


(* ::Text::Initialization:: *)
(*(*We will then be able to generate all our desired coefficients for all Cutksoky cuts as follows:*)*)


Print["Now starting the computation of coefficients..."]


(* ::Input::Initialization:: *)
ComputeCoefficientsForCutkoskyCut[iCutkosky_,allCoefficientsForThisCutkoskyNum_,chosenTermIndexOffset_]:=Block[{tmp,FinalAllCoefficients},
tmp=PrintTemporary["Now generating coefficients for cutkosky cut #"<>ToString[iCutkosky]<>" ..."];
FinalAllCoefficients=<|
Keys[FinalAllSelectedCutkoskyCutsMomentumBases][[iCutkosky]]->GenerateNumeratorCoefficients[
allCoefficientsForThisCutkoskyNum,
DEBUG->True,
msgPrefix->("Kernel ID="<>ToString[$KernelID]<>" Cut #"<>ToString[iCutkosky]<>" : "),
termIndexOffset->chosenTermIndexOffset,
printoutFrequency->10,
ReplacementRules->{},
consistencyCheckLevel->0
]
|>;
NotebookDelete[tmp];
FinalAllCoefficients
]


(* ::Text:: *)
(*Utilities for performing the computation*)


(* ::Input::Initialization:: *)
chunks[list_,nChunks_,chunkID_]:=Block[{chunkSize=QuotientRemainder[Length[list],nChunks][[1]]+1},
{
(chunkID-1) chunkSize+1,
Min[chunkID chunkSize,Length[list]]
}
]


(* ::Input::Initialization:: *)
GetAllCoeffs[iCutkosky_,MinTermIndex_,MaxTermIndex_,OptionsPattern[{Parallel->False}]]:=Block[
{NKernels,allCoefficientsForThisCutkoskyNum,allParallelCoefsSummed={},allParallelEvaluations,i},

allCoefficientsForThisCutkoskyNum=List@@(Expand[(ExpandScalarProduct[
Total[ProcessedSelectedGraphNumeratorAsList[[MinTermIndex;;MaxTermIndex]]]/.FinalAllSelectedCutkoskyCutsMomentumBases[[iCutkosky]]
])/.{Pair[Momentum[p1],Momentum[p2]]:>2 10^6}]);

(*Sequential evaluation*)
If[Not[OptionValue[Parallel]],Return[ComputeCoefficientsForCutkoskyCut[iCutkosky,allCoefficientsForThisCutkoskyNum,0]]];

(* Parallel evaluation below *)
NKernels=Max[ParallelEvaluate[$KernelID]];
(* Load context in Parallel kernels *)
ParallelEvaluate[
SetDirectory[NoteBookDir];
Get[LTDToolsPath];
];

allParallelEvaluations=ParallelEvaluate[
Block[{chunk,printoutIndexOffset},
chunk=chunks[allCoefficientsForThisCutkoskyNum,NKernels,$KernelID];
printoutIndexOffset=(chunk[[1]]-1);
ComputeCoefficientsForCutkoskyCut[iCutkosky,allCoefficientsForThisCutkoskyNum[[chunk[[1]];;chunk[[2]]]],printoutIndexOffset]
]
];
<|Keys[FinalAllSelectedCutkoskyCutsMomentumBases][[iCutkosky]]->SumCoeffs[Table[c[[1]],{c,allParallelEvaluations}]]|>
];


(* ::Text:: *)
(*Now finally perform the computation*)


(* ::Input::Initialization:: *)
FinalAllCoefficients=<||>;
Block[{iCutkosky,tmp,
(* Select below the range of terms to consider *)
MinTermIndex=1,
MaxTermIndex=-1,
MaxiCutksoksky=Length[FinalAllSelectedCutkoskyCutsMomentumBases]
},
For[iCutkosky=1,iCutkosky<=MaxiCutksoksky,iCutkosky++,
AppendTo[FinalAllCoefficients,GetAllCoeffs[iCutkosky,MinTermIndex,MaxTermIndex,Parallel->True]];
Export["AllCoefficients_up_to_cut_"<>ToString[iCutkosky]<>".m",FinalAllCoefficients];
];
];
Export["FinalAllCoefficients.m",FinalAllCoefficients]
