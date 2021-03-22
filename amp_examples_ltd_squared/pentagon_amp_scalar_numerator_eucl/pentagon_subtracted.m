(* ::Package:: *)

Quit[]


andrea = {
a3->-((p1.q1+p1.q2)^2/(p1.p2 (m1+m2-2 p1.q1-2 p1.q2+2 q1.q2) (m1 p1.q2-p1.q1 (m2+2 q1.q2)))),
a4->(p1.q1)^2/(p1.p2 (m1-2 p1.q1) (m1 p1.q2-p1.q1 (m2+2 q1.q2))),
a8->(-p1.p2+p2.q1+p2.q2)^2/(p1.p2 (m1+m2-2 p1.q1-2 p1.q2+2 q1.q2) ((m1-2 p1.q1) p2.q2+p1.p2 (m2-2 p1.q2+2 q1.q2)-p2.q1 (m2-2 p1.q2+2 q1.q2))),
a9->-((p1.p2-p2.q1)^2/(p1.p2 (m1-2 p1.q1) ((m1-2 p1.q1) p2.q2+p1.p2 (m2-2 p1.q2+2 q1.q2)-p2.q1 (m2-2 p1.q2+2 q1.q2)))),
a10->1/((m1-2 p1.q1) (m1+m2-2 p1.q1-2 p1.q2+2 q1.q2))
};






rr = {p1 -> (Rationalize[#,0]&/@{1, 0, 0, 1}),
p2 -> (Rationalize[#,0]&/@{-1, 0., 0., 1.}),
q1 -> (Rationalize[#,0]&/@{0.,2.3,0.2,1.1}),
q2 -> (Rationalize[#,0]&/@{0.,0.2,1.3,3.123})};


andreaNum=(andrea /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3]) /. rr /. sp[a_,b_]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]] //N


andreaNumEucl=(andrea /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3]) /. rr2 /. sp[a_,b_]:>a[[1]]*b[[1]] - Dot[a[[2;;]],b[[2;;]]] //N


rrZenoToAndrea =
Sort@(Rule[-Values[#],-Keys[#]]&/@{ a3 -> -xZenoConst4,
 a4 -> -xZenoConst5,
 - a8 -> - xZenoConst1,
 a9 -> -xZenoConst2,
 a10 -> -xZenoConst3
})


rrZenoToAndrea /. andreaNumEucl


numerator = 1-xZenoConst1*(k1.k1+2*k1.p1+p1.p1)*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)+xZenoConst2*(k1.k1+2*k1.p1+p1.p1)*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)+xZenoConst3*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)+xZenoConst4*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)*(k1.k1-2*k1.p2+p2.p2)+xZenoConst5*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)*(k1.k1-2*k1.p2+p2.p2);


StringReplace[#,{"["->"(","]"->")"}]&@(ToString[(numerator /. Dot -> sp),InputForm])


test0 = numerator /. setTo0/. rrZenoToAndrea /. andrea  /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3] /. q1->p3 /. q2->p4 /. q3 ->p5 /. k3 ->p4 /. k2 -> p3//Simplify;


setTo0 = {xZenoConst1 ->0, xZenoConst2 ->0, xZenoConst4->0,xZenoConst5 ->0};


test1 = StringReplace[#,{"["->"(","]"->")"}]&@(ToString[#,InputForm]&@(numerator-1 /. setTo0 /. rrZenoToAndrea /. andrea /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3] /. q1->p3 /. q2->p4 /. q3 ->p5 /. k3 ->p4 /. k2 -> p3//Simplify)) ;


test1b = StringReplace[#,{"["->"(","]"->")"}]&@(ToString[#,InputForm]&@(numerator-1 /. setTo0   /. Dot->sp /. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3] /. q1->p3 /. q2->p4 /. q3 ->p5 /. k3 ->p4 /. k2 -> p3//Simplify)) ;


test1b


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


test1


(1-test0 /. sp -> Dot /. p3 -> k2 /. p4 -> k3)/ampDenomTest//FullSimplify


test2//Simplify


ampDenomTest


a10 /.andrea


rrZenoToAndrea /. andrea


rrZenoToAndrea /. andrea/. Dot -> sp/. m1 ->sp[q1,q1] /. m2 -> sp[q2,q2] /. m3 -> sp[q3,q3] /. q1->k2 /. q2->k3 /. q3 ->p5 /. sp -> Dot


ampDenomTest
