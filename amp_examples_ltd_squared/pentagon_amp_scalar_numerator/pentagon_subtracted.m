(* ::Package:: *)

Quit[]


andrea = {
a3->-((p1.q1+p1.q2)^2/(p1.p2 (m1+m2-2 p1.q1-2 p1.q2+2 q1.q2) (m1 p1.q2-p1.q1 (m2+2 q1.q2)))),
a4->(p1.q1)^2/(p1.p2 (m1-2 p1.q1) (m1 p1.q2-p1.q1 (m2+2 q1.q2))),
a8->(-p1.p2+p2.q1+p2.q2)^2/(p1.p2 (m1+m2-2 p1.q1-2 p1.q2+2 q1.q2) ((m1-2 p1.q1) p2.q2+p1.p2 (m2-2 p1.q2+2 q1.q2)-p2.q1 (m2-2 p1.q2+2 q1.q2))),
a9->-((p1.p2-p2.q1)^2/(p1.p2 (m1-2 p1.q1) ((m1-2 p1.q1) p2.q2+p1.p2 (m2-2 p1.q2+2 q1.q2)-p2.q1 (m2-2 p1.q2+2 q1.q2)))),
a10->1/((m1-2 p1.q1) (m1+m2-2 p1.q1-2 p1.q2+2 q1.q2))
};



rr={p1 -> (Rationalize[#,0]&/@{500, 0, 0, 500}),
p2 -> (Rationalize[#,0]&/@{500., 0., 0., -500.}),
q1 -> (Rationalize[#,0]&/@{349.76387887, 197.3970875, 277.2935156, -4.06480273}),
q2 -> (Rationalize[#,0]&/@{274.3666484, 94.14456732, -156.29256755, -183.49625453}),
q3 -> (Rationalize[#,0]&/@{375.86947272, -291.54165482, -121.00094805, 187.56105726})};


andreaNum=(andrea /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3]) /. rr /. sp[a_,b_]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]];


rrZenoToAndrea =
Sort@(Rule[-Values[#],-Keys[#]]&/@{ a3 -> -xZenoConst4,
 a4 -> -xZenoConst5,
 - a8 -> - xZenoConst1,
 a9 -> -xZenoConst2,
 a10 -> -xZenoConst3
})


numerator = 1-xZenoConst1*(k1.k1+2*k1.p1+p1.p1)*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)+xZenoConst2*(k1.k1+2*k1.p1+p1.p1)*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)+xZenoConst3*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)+xZenoConst4*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)*(k1.k1-2*k1.p2+p2.p2)+xZenoConst5*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)*(k1.k1-2*k1.p2+p2.p2);


test0 = numerator /. rrZenoToAndrea /. andrea  /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3] /. q1->p3 /. q2->p4 /. q3 ->p5 /. k3 ->p4 /. k2 -> p3//Simplify;


(* this is the input for pentagon_new *)
test1 = StringReplace[#,{"["->"(","]"->")"}]&@(ToString[#,InputForm]&@(numerator /. rrZenoToAndrea /. andrea  /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3] /. q1->p3 /. q2->p4 /. q3 ->p5 /. k3 ->p4 /. k2 -> p3//Simplify)) ;


test1


(* FROM HERE ONWARDS ITS JUST CHECKING STUFF *)
(* after wrapping everything into the denominator wrapper with FORM *)
test2 =(ToExpression[#,TraditionalForm]&/@(ToExpression[StringReplace[Import["./TEST_pentagon_new/process_dir/FORM/workspace/SG_0_color_decomp.txt"],{"\"color(1)\":"->""}],TraditionalForm]))[[1]];


(* check the content of the denom wrapper *)
ampDenomTest =(Cases[test2,_ampDenom,Infinity] //Union) /. Dot ->sp /. sp[a_,b_^c_]:>sp[a,b]^c /. sp :>Dot /. ampDenom[x_]:>1/x //FullSimplify;


andreaTransformed = 1/(Denominator/@(Values/@(andrea /. m1 -> k2.k2 /. m2 -> k3.k3 /. q1 -> k2 /. q2->k3)))//FullSimplify;
andreaTransformed2 = ((Values/@(andrea /. m1 -> k2.k2 /. m2 -> k3.k3 /. q1 -> k2 /. q2->k3)))//FullSimplify;


ampDenomTest/. Dot ->sp /. sp[a_,b_^c_]:>sp[a,b]^c /. ampDenom[x_]:>1/x/. k2 ->q1 /. k3 -> q2 /. p3->q1 /. p4->q2 /. p5 -> q3 /. rr /. sp[a_List,b_List]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]] // FullSimplify //N //Sort
andreaTransformed /. Dot ->sp /. sp[a_,b_^c_]:>sp[a,b]^c /. ampDenom[x_]:>1/x/. k2 ->q1 /. k3 -> q2 /. p3->q1 /. p4->q2 /. p5 -> q3 /. rr /. sp[a_List,b_List]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]] // FullSimplify //N //Sort


andreaTransformed2 /. Dot ->sp /. sp[a_,b_^c_]:>sp[a,b]^c /. ampDenom[x_]:>1/x/. k2 ->q1 /. k3 -> q2 /. p3->q1 /. p4->q2 /. p5 -> q3 /. rr /. sp[a_List,b_List]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]] // FullSimplify //N //Sort
Values@andreaNum //N //Sort


(* check the complete expression *)
tt2 = test2/. Dot ->sp /. sp[a_,b_^c_]:>sp[a,b]^c /. ampDenom[x_]:>1/x /. k2 ->q1 /. k3 -> q2  /. k1->{kk1,kk2,kk3,kk4} /. rr /. sp[a_List,b_List]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]] //Expand;
tt2-ToExpression[test1 ,TraditionalForm]/. k2 ->q1 /. k3 -> q2 /. p3->q1 /. p4->q2 /. p5 -> q3 /. k1->{kk1,kk2,kk3,kk4} /. rr /. sp[a_List,b_List]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]] // FullSimplify
