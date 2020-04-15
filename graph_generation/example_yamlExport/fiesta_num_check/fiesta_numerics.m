(* ::Package:: *)

(* ::Section:: *)
(*(*Compute singular behaviour with FIESTA *)*)


(* ::Input::Initialization:: *)

FIESTAPath="/home/arminsch/my_programs/fiesta/FIESTA4/";
Get[FIESTAPath<>"/FIESTA4.m"]



NumberOfSubkernels=36;
NumberOfLinks=36;
MPThreshhold=100;
ReturnErrorWithBrackets=True;
CurrentIntegratorSettings={{"epsrel","1.000000E-08"},{"epsabs","1.000000E-12"},{"mineval","1000000"},{"maxeval","5000000"},{"nstart","100000"},{"nincrease","500"},{"seed","0"},{"rng","0"}};
DigitsLimit=10;
UsingC=True;

(* Load Master *)
kinematics=<|"mphi"->125,
Thread@Rule[{p1,p2,p3,p4},{1,1,-1,-1}*SetPrecision[ImportString["[[6.0,0.0,0.0,5.91607978309962],
[6.0,0.0,0.0,-5.91607978309962],
[-6.0,-1.3124738333059,-5.26330888118183,2.36114210884473],
[-6.0,1.3124738333059,5.26330888118183,-2.36114210884473]]","PythonExpression"],32]]
|>;
props={{-k1^2+mphi^2,mphi^2-(k1+p3)^2,-k2^2+mphi^2,mphi^2-(k2+p4)^2,-(k1+k2)^2+mphi^2,mphi^2-(-k1-k2-p1-p2)^2},{-k1^2+mphi^2,mphi^2-(k1-p2)^2,-k2^2+mphi^2,mphi^2-(k2+p4)^2,-(k1+k2)^2+mphi^2,mphi^2-(-k1-k2-p1+p3)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,-(k1+k2)^2+mphi^2,mphi^2-(-k1-k2+p3+p4)^2},{-k1^2+mphi^2,mphi^2-(k1+p3)^2,-k2^2+mphi^2,mphi^2-(k2+p4)^2,mphi^2-(-k1-p1-p2)^2,mphi^2-(k1-k2+p3)^2},{-k1^2+mphi^2,mphi^2-(k1-p2)^2,-k2^2+mphi^2,mphi^2-(k2+p4)^2,mphi^2-(-k1-p1+p3)^2,mphi^2-(k1-k2-p2)^2},{-k1^2+mphi^2,mphi^2-(k1+p4)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,mphi^2-(-k1-p1+p3)^2,mphi^2-(k1-k2+p4)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,mphi^2-(-k1+p3+p4)^2,mphi^2-(k1-k2-p1)^2},{-k1^2+mphi^2,mphi^2-(k1+p3)^2,mphi^2-(-k1+p4)^2,-k2^2+mphi^2,mphi^2-(k2-p1-p2)^2,mphi^2-(k1-k2+p3)^2},{-k1^2+mphi^2,mphi^2-(k1-p2)^2,mphi^2-(-k1+p4)^2,-k2^2+mphi^2,mphi^2-(k2-p1+p3)^2,mphi^2-(k1-k2-p2)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1-p2)^2,-k2^2+mphi^2,mphi^2-(k2+p3+p4)^2,mphi^2-(k1-k2-p1)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1-p2)^2,-k2^2+mphi^2,mphi^2-(k2+p3)^2,mphi^2-(-k1+k2+p1)^2,mphi^2-(k1-k2-p1+p4)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1+p3)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,mphi^2-(-k1+k2+p1)^2,mphi^2-(k1-k2-p1+p4)^2},{-k1^2+mphi^2,mphi^2-(k1+p3)^2,mphi^2-(-k1+p4)^2,-k2^2+mphi^2,mphi^2-(k2-p1)^2,mphi^2-(-k1+k2-p3)^2,mphi^2-(k1-k2-p2+p3)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1-p2)^2,-k2^2+mphi^2,mphi^2-(k2+p3)^2,mphi^2-(-k2+p4)^2,mphi^2-(k1+k2-p1+p3)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1+p3)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,mphi^2-(-k2+p4)^2,mphi^2-(k1+k2-p1-p2)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1+p3)^2,-k2^2+mphi^2,mphi^2-(k2+p4)^2,mphi^2-(-k2-p2)^2,mphi^2-(k1+k2-p1+p4)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1-p2)^2,-k2^2+mphi^2,mphi^2-(k2+p3)^2,mphi^2-(-k2+p4)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1+p3)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,mphi^2-(-k2+p4)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1-p2)^2,mphi^2-(k1-p1+p3)^2,-k2^2+mphi^2,mphi^2-(k2+p4)^2,mphi^2-(-k1-k2-p2)^2},{-k1^2+mphi^2,mphi^2-(k1-p1)^2,mphi^2-(-k1+p3)^2,mphi^2-(k1-p1+p4)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,mphi^2-(-k1-k2+p3)^2},{-k1^2+mphi^2,mphi^2-(k1+p3)^2,mphi^2-(-k1-p1)^2,mphi^2-(k1-p2+p3)^2,-k2^2+mphi^2,mphi^2-(k2+p4)^2,mphi^2-(-k1-k2-p1)^2},{-k1^2+mphi^2,mphi^2-(k1+p3)^2,mphi^2-(-k1-p1)^2,mphi^2-(k1+p3+p4)^2,-k2^2+mphi^2,mphi^2-(k2-p2)^2,mphi^2-(-k1-k2-p1)^2}};
SetAttributes[sp,Orderless]
sp[x_List,y_List]:=x[[1]]y[[1]]-x[[2;;]].y[[2;;]]


replRule=Join[Thread@Rule[Times@@@(Subsets[{p1,p2,p3,p4,p1,p2,p3,p4},{2}]//Union),(sp[#[[1]],#[[2]]]&/@(Subsets[{p1,p2,p3,p4,p1,p2,p3,p4},{2}]//Union)/.kinematics)],kinematics[[{"mphi"}]]//Normal];


(* uncomment if you rerun is neccessary *)
Do[
count=count+1;
ReturnErrorWithBrackets=True;
ComplexMode=True;
Print["diag "<>ToString@count<>" started at"<>DateString[]]
int1=(-1)^Length@diag*SDEvaluate[UF[{k1,k2},diag,(replRule)],ConstantArray[1,Length@diag],0];
Export["./ps1/diag"<>ToString[count]<>".m",{int1}];
Print["diag "<>ToString@count<>" finished at"<>DateString[]]
,{diag,props}]
 

