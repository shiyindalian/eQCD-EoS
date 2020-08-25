(* ::Package:: *)

Clear["`*"]


lF00=1/3 (1-etapsi/4) 1/Sqrt[1+mf2[rho,rho2]] (1-nfa[k Sqrt[1+mf2[rho,rho2]]]-nff[k Sqrt[1+mf2[rho,rho2]]]);


lF10=D[lF00,rho];
lF20=D[lF10,rho];
lF30=D[lF20,rho];
lF40=D[lF30,rho];
lF50=D[lF40,rho];
lF01=D[lF00,rho2];
lF11=D[lF01,rho];
lF21=D[lF11,rho];
lF31=D[lF21,rho];
lF02=D[lF01,rho2];
lF12=D[lF02,rho];


repla={nfa[k Sqrt[1+mf2[rho,rho2]]]->nfad0x,
\!\(\*SuperscriptBox[\(nfa\), 
TagBox[
RowBox[{"(", "1", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nfad1x,
\!\(\*SuperscriptBox[\(nfa\), 
TagBox[
RowBox[{"(", "2", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nfad2x,
\!\(\*SuperscriptBox[\(nfa\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nfad3x,
\!\(\*SuperscriptBox[\(nfa\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nfad4x,
\!\(\*SuperscriptBox[\(nfa\), 
TagBox[
RowBox[{"(", "5", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nfad5x,nff[k Sqrt[1+mf2[rho,rho2]]]->nffd0x,
\!\(\*SuperscriptBox[\(nff\), 
TagBox[
RowBox[{"(", "1", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nffd1x,
\!\(\*SuperscriptBox[\(nff\), 
TagBox[
RowBox[{"(", "2", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nffd2x,
\!\(\*SuperscriptBox[\(nff\), 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nffd3x,
\!\(\*SuperscriptBox[\(nff\), 
TagBox[
RowBox[{"(", "4", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nffd4x,
\!\(\*SuperscriptBox[\(nff\), 
TagBox[
RowBox[{"(", "5", ")"}],
Derivative],
MultilineFunction->None]\)[k Sqrt[1+mf2[rho,rho2]]]->nffd5x,
mf2[rho,rho2]->mf2d00rho,
\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mf2d10rho,
\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mf2d01rho,
\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]->mf2d02rho};



\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"3", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"4", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"5", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"3", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"4", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"3", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "4"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "4"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;

\!\(\*SuperscriptBox[\(mf2\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "5"}], ")"}],
Derivative],
MultilineFunction->None]\)[rho,rho2]=0;


lF00=lF00/.repla;
lF10=lF10/.repla;
lF20=lF20/.repla;
lF30=lF30/.repla;
lF40=lF40/.repla;
lF50=lF50/.repla;
lF01=lF01/.repla;
lF11=lF11/.repla;
lF21=lF21/.repla;
lF31=lF31/.repla;
lF02=lF02/.repla;
lF12=lF12/.repla;


FortranForm[{
lft00==lF00,
lft10==lF10,
lft20==lF20,
lft30==lF30,
lft40==lF40,
lft50==lF50,
lft01==lF01,
lft11==lF11,
lft21==lF21,
lft31==lF31,
lft02==lF02,
lft12==lF12
}]>>"./lFm.f90";
