(* ::Package:: *)

Clear["`*"]


lB00=2/3 (1-etaphi/ 5) 1/ Sqrt[1+mb2[rho,rho2]] (1/ 2+nbo[k Sqrt[1+mb2[rho,rho2]]]);


lB10=D[lB00,rho];
lB20=D[lB10,rho];
lB30=D[lB20,rho];
lB40=D[lB30,rho];
lB50=D[lB40,rho];
lB01=D[lB00,rho2];
lB11=D[lB01,rho];
lB21=D[lB11,rho];
lB31=D[lB21,rho];
lB02=D[lB01,rho2];
lB12=D[lB02,rho];


repla={nbo[k Sqrt[1+mb2[rho,rho2]]]->nbd0x,
\!\(\*SuperscriptBox[\(nbo\), 
TagBox[
RowBox[{"(", "1", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mb2[rho,rho2]]]->nbd1x,
\!\(\*SuperscriptBox[\(nbo\), 
TagBox[
RowBox[{"(", "2", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mb2[rho,rho2]]]->nbd2x,
\!\(\*SuperscriptBox[\(nbo\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mb2[rho,rho2]]]->nbd3x,
\!\(\*SuperscriptBox[\(nbo\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mb2[rho,rho2]]]->nbd4x,
\!\(\*SuperscriptBox[\(nbo\), 
TagBox[
RowBox[{"(", "5", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mb2[rho,rho2]]]->nbd5x,mb2[rho,rho2]->mb2d00rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d10rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d20rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"3", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d30rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"4", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d40rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"5", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d50rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d01rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d11rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d21rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"3", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d31rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d02rho,
\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mb2d12rho};



\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"4", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"3", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "4"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "4"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mb2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "5"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;


lB00=lB00/.repla;
lB10=lB10/.repla;
lB20=lB20/.repla;
lB30=lB30/.repla;
lB40=lB40/.repla;
lB50=lB50/.repla;
lB01=lB01/.repla;
lB11=lB11/.repla;
lB21=lB21/.repla;
lB31=lB31/.repla;
lB02=lB02/.repla;
lB12=lB12/.repla;


FortranForm[{
lbt00==lB00,
lbt10==lB10,
lbt20==lB20,
lbt30==lB30,
lbt40==lB40,
lbt50==lB50,
lbt01==lB01,
lbt11==lB11,
lbt21==lB21,
lbt31==lB31,
lbt02==lB02,
lbt12==lB12
}]>>"lbm.f90";
