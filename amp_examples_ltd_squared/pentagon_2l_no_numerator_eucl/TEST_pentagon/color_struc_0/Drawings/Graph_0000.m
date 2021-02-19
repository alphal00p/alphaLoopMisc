GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=Labeled[GraphClass[{
Labeled[Style[CreateEdge["1","3",1],Gray,Thickness[0.002000]],"psis|q1|p1"],
Labeled[Style[CreateEdge["2","4",2],Gray,Thickness[0.002000]],"psis|q2|p2"],
Labeled[Style[CreateEdge["10","8",14],Gray,Thickness[0.002000]],"psis|q3|p1"],
Labeled[Style[CreateEdge["10","9",15],Gray,Thickness[0.002000]],"psis|q4|p2"],
Labeled[Style[CreateEdge["4","11",3],Gray,Thickness[0.005000]],"psis|pq1|k1"],
Labeled[Style[CreateEdge["3","11",4],Gray,Thickness[0.005000]],"psis|pq2|k2"],
Labeled[Style[CreateEdge["11","3",5],Gray,Thickness[0.002000]],"psis|pq3|k1+k2"],
Labeled[Style[CreateEdge["3","5",6],Gray,Thickness[0.002000]],"psis|pq4|k1+p1"],
Labeled[Style[CreateEdge["5","6",7],Gray,Thickness[0.002000]],"psis|pq5|k1-k3+p1"],
Labeled[Style[CreateEdge["6","7",8],Gray,Thickness[0.002000]],"psis|pq6|k1-k3-k4+p1"],
Labeled[Style[CreateEdge["7","4",9],Gray,Thickness[0.002000]],"psis|pq7|k1-p2"],
Labeled[Style[CreateEdge["5","10",10],Gray,Thickness[0.005000]],"psis|pq8|k3"],
Labeled[Style[CreateEdge["6","10",11],Gray,Thickness[0.005000]],"psis|pq9|k4"],
Labeled[Style[CreateEdge["7","10",12],Gray,Thickness[0.002000]],"psis|pq10|-k3-k4+p1+p2"]
},
EdgeShapeFunction -> {
CreateEdge["1","3",1]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["2","4",2]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["10","8",14]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["10","9",15]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["4","11",3]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["3","11",4]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["11","3",5]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["3","5",6]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["5","6",7]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["6","7",8]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["7","4",9]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["5","10",10]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["6","10",11]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["7","10",12]->GraphElementData["Arrow", "ArrowSize" -> 0.015000]
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 8,Small],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 8, Small],
VertexSize -> small,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> {"SpringEmbedding"},
ImageSize -> {660.000000, 510.000000}
],"MG: SG_QG0 | FORM: #0"];
Export["Graph_0000.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {825.000000, 637.500000}];