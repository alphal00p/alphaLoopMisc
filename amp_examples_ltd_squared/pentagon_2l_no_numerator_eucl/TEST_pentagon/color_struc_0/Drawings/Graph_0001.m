GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=Labeled[GraphClass[{
Labeled[Style[CreateEdge["1","3",1],Gray,Thickness[0.002000]],"psis|q1|p1"],
Labeled[Style[CreateEdge["2","4",2],Gray,Thickness[0.002000]],"psis|q2|p2"],
Labeled[Style[CreateEdge["10","8",13],Gray,Thickness[0.002000]],"psis|q3|p1"],
Labeled[Style[CreateEdge["10","9",14],Gray,Thickness[0.002000]],"psis|q4|p2"],
Labeled[Style[CreateEdge["4","3",3],Gray,Thickness[0.005000]],"psis|pq1|k1"],
Labeled[Style[CreateEdge["3","3",4],Gray,Thickness[0.005000]],"psiuv|pq2|k2"],
Labeled[Style[CreateEdge["3","5",5],Gray,Thickness[0.002000]],"psis|pq3|k1+p1"],
Labeled[Style[CreateEdge["5","6",6],Gray,Thickness[0.002000]],"psis|pq4|k1-k3+p1"],
Labeled[Style[CreateEdge["6","7",7],Gray,Thickness[0.002000]],"psis|pq5|k1-k3-k4+p1"],
Labeled[Style[CreateEdge["7","4",8],Gray,Thickness[0.002000]],"psis|pq6|k1-p2"],
Labeled[Style[CreateEdge["5","10",9],Gray,Thickness[0.005000]],"psis|pq7|k3"],
Labeled[Style[CreateEdge["6","10",10],Gray,Thickness[0.005000]],"psis|pq8|k4"],
Labeled[Style[CreateEdge["7","10",11],Gray,Thickness[0.002000]],"psis|pq9|-k3-k4+p1+p2"]
},
EdgeShapeFunction -> {
CreateEdge["1","3",1]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["2","4",2]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["10","8",13]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["10","9",14]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["4","3",3]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["3","3",4]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["3","5",5]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["5","6",6]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["6","7",7]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["7","4",8]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["5","10",9]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["6","10",10]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["7","10",11]->GraphElementData["Arrow", "ArrowSize" -> 0.015000]
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 8,Small],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 8, Small],
VertexSize -> small,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> {"SpringEmbedding"},
ImageSize -> {660.000000, 510.000000}
],"MG: SG_QG1 | FORM: #1"];
Export["Graph_0001.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {825.000000, 637.500000}];