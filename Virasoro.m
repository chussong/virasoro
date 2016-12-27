(* ::Package:: *)

BeginPackage["Virasoro`"]

(*For simple computations*)
h::usage = " h[m_,n_] gives degenerate Virasoro operator dimension h_mn.";
\[Lambda]sq::usage = "\[Lambda]sq[p_,q_] gives square of internal \[Lambda]_pq.";
\[Lambda]::usage = "\[Lambda][p_,q_] gives internal \[Lambda]_pq.";
Bsq::usage = "Bsq[c_] gives b^2 corresponding to the given c.";
c::usage = "c[bsq_] gives c corresponding to the given b^2.";

(*For plotting in Lorentzian time*)
EKMono::usage = "EKMono[r_,tL_] gives the EllipticK[z] for these parameters accounting for the monodromy.";
qVal::usage = "qVal[r_,tL] gives the q value corresponding to the given Lorentzian time tL using radius r.";
VofT::usage = "QtoT[coeVec_,c_,hL_,hH_,hp_,r_,tL_,maxOrder_] gives the Lorentzian Virasoro block approximated by coeVec at Lorentzian time tL.";
SemiClassical::usage = "SemiClassical[c_,hL_,hH_,r_,tL_] gives the value of the semiclassical approximation of the vacuum block at Lorentzian time tL.";
DegenBlock12::usage = "DegenBlock12[c_,hL_,hh_,r_,tL_] gives the exact degenerate h_12 vacuum block at time tL.";
DegenBlock21::usage = "DegenBlock21[c_,hL_,hh_,r_,tL_] gives the exact degenerate h_21 vacuum block at time tL.";

(*For reading results from virasoro batch runs*)
VLengthCheck::usage = "VLengthCheck[filename_] gives the length of the result vector filename.";
VRead::usage = "VRead[filename_] gives the run results as a vector of vectors of numbers, ready to be used in Mathematica.";
VPlot::usage = "VPlot[results_,startingRun_,runStep_,r_(,Option\[Rule]Value)] plots all runs in the VRead output results, using radius r. runStep>1 allows you to skip runs.

Options are specified as in built-in Plot: Option\[Rule]Value, e.g. EndTime\[Rule]40. Options are EndTime, PointsPerTL, and Compare. Setting Compare\[Rule]\"Semi\" compares to semiclassical vacuum block; other possibilities are \"12\" and \"21\" which compare to exact degenerate blocks.";

Begin["Private`"]
h[m_,n_]:=b^2*(1-n^2)/4+(1-m^2)/(4*b^2)+(1-m*n)/2;
\[Lambda]sq[p_,q_]:=4*h[p,q]-(b+1/b)^2;
\[Lambda][p_,q_]:=Sqrt[\[Lambda]sq[p,q]];
Bsq[c_]:=(-13+c+Sqrt[25-26c+c^2])/12;
c[bsq_]:=13+6*(bsq+1/bsq);

EKMono[r_,tL_]:=EllipticK[1-r*Exp[-I*tL]]-2*I*(1+Floor[(-tL-\[Pi])/(2\[Pi])]) EllipticK[r*Exp[-I*tL]];
qVal[r_,tL_]:=Exp[-\[Pi]*EllipticK[r*Exp[-I*tL]]/EKMono[r,tL]];
VofT[coeVec_,c_,hL_,hH_,hp_,r_,tL_,maxOrder_]:=Module[{z,q},
q=qVal[r,tL];
z=1-r*E^(-I tL);
(r^hL E^(-I hL (tL+3\[Pi])))(16 q)^(hp-(c-1)/24) z^((c-1)/24-2hL) (r)^((c-1)/24-hH-hL) E^(-I tL ((c-1)/24-hH-hL)) EllipticTheta[3,0,q]^((c-1)/2-8(hH+hL))*(Table[q^i,{i,0,maxOrder,2}].Take[coeVec,1+Floor[maxOrder/2]])
];
SemiClassical[c_,hL_, hH_,r_,tL_]:=Module[{\[Alpha]},(r^hL E^(-I hL (tL+3\[Pi])))((\[Alpha]^(2 hL) r^((\[Alpha]-1)hL) E^(-I tL(\[Alpha]-1)hL))/(1-r^\[Alpha] E^(-I tL \[Alpha]))^(2hL))/.\[Alpha]->Sqrt[1-(24hH)/c]
];
DegenBlock12[c_,hL_,hh_,r_,tL_]:=Module[{a1,b1,c1,z,bsq},
bsq=1/Bsq[c];
1/(1-r*Exp[-I tL])^(2hL) r^b1 Exp[-I b1 tL/2]((Gamma[a1+b1-c1] Gamma[c1])/(Gamma[a1] Gamma[b1]) r^(-a1-b1+c1) Exp[-I (-tL a1-tL b1+tL c1)] Hypergeometric2F1[-a1+c1,-b1+c1,1-a1-b1+c1,r E^(-I tL)] +(Gamma[-a1-b1+c1] Gamma[c1])/(Gamma[c1-a1] Gamma[c1-b1]) Hypergeometric2F1[a1,b1,a1+b1-c1+1,r E^(-I tL)] )/.{a1->1+1/bsq,b1->(1+bsq+Sqrt[1+bsq^2+bsq *(2-4 hh)])/bsq,c1->2+2/bsq}/.z->1-r E^(-I tL)
];
DegenBlock21[c_,hL_,hh_,r_,tL_]:=Module[{a1,b1,c1,z,bsq},
bsq=Bsq[c];
1/(1-r Exp[-I tL])^(2hL) r^b1 Exp[-I b1 tL/2]((Gamma[a1+b1-c1] Gamma[c1])/(Gamma[a1] Gamma[b1]) r^(-a1-b1+c1) Exp[-I (-tL a1-tL b1+tL c1)] Hypergeometric2F1[-a1+c1,-b1+c1,1-a1-b1+c1,r E^(-I tL)] +(Gamma[-a1-b1+c1] Gamma[c1])/(Gamma[c1-a1] Gamma[c1-b1]) Hypergeometric2F1[a1,b1,a1+b1-c1+1,r E^(-I tL)] )/.{a1->1+bsq,b1->(1+1/bsq+Sqrt[1+bsq^(-2)+1/bsq *(2-4 hh)])*bsq,c1->2+2*bsq}/.z->1-r E^(-I tL)
];

VLengthCheck[filename_]:=Module[{resultFile,entries, null},
resultFile=OpenRead[filename];
entries =0;
null[_String]:=Null;
entries=Length@ReadList[filename,null@String,NullRecords->True];
Close[filename];
Return[entries];
];
VRead[filename_]:=Module[{resultFile,entries,null,results},
entries=VLengthCheck[filename];
results=ConstantArray[0,entries];
resultFile=OpenRead[filename];
Do[
results[[2i-1]]=ToExpression@ReadLine[resultFile];
Do[results[[2i-1]][[j]]=ToString[results[[2i-1]][[j]]],{j,1,5}];
results[[2i]]=ToExpression@ReadLine[resultFile];
,{i,1,entries/2}];
Close[filename];
Return[results];
];
Options[VPlot]={EndTime->40, Compare->{}, PointsPerTL->1};
VPlot[results_,startingRun_,runStep_,r_,OptionsPattern[]]:=Module[{c,hl,hh,hp,endTime, compareVec, plotVector, plotLegends},
endTime=OptionValue[EndTime];
If[VectorQ[OptionValue[Compare]],compareVec = OptionValue[Compare], compareVec = {OptionValue[Compare]}];
Do[If[Length@results[[2*i]]<10,Continue[]];
c=results[[2i-1]][[1]];
hl=results[[2i-1]][[2]];
hh=results[[2i-1]][[3]];
hp=results[[2i-1]][[4]];
Print[ListLogLogPlot[{results[[2*i]],-results[[2*i]]},PlotLabel->"c="<>StringTake[c,Min[5,StringLength[c]]]<>"    \!\(\*SubscriptBox[\(h\), \(l\)]\)="<>StringTake[hl,Min[5,StringLength[hl]]]<>"    \!\(\*SubscriptBox[\(h\), \(h\)]\)="<>StringTake[hh,Min[5,StringLength[hh]]]<>"    \!\(\*SubscriptBox[\(h\), \(p\)]\)="<>StringTake[hp,Min[5,StringLength[hp]]],PlotMarkers->".",PlotStyle->{Lighter@Blue,Lighter@Red},PlotLegends->{"Positive","Negative"},PlotLabel->"Coefficient values in q space",DataRange->{0,2*Length@results[[2*i]]}]];
(*Print[ListPlot[LogFluct[results[[2*i]]],PlotLabel\[Rule]"Fluctuations about smoothed average log", PlotMarkers\[Rule]".",ColorFunction\[Rule]Coloring,ColorFunctionScaling\[Rule]False,DataRange\[Rule]{0,2*Length@results[[2*i]]}]];*)
plotVector := {Abs@VofT[results[[2*i]],ToExpression@c,ToExpression@hl,ToExpression@hh,ToExpression@hp,r,tL,2*Length@results[[2*i]]-2]};
plotLegends := {"Computed"};
If[MemberQ[compareVec, "Semi"],
plotVector=Append[Abs@SemiClassical[ToExpression@c,ToExpression@hl,ToExpression@hh,r,tL]]@plotVector;
plotLegends=Append["Semiclassical"]@plotLegends;
];
If[MemberQ[compareVec, "12"],
plotVector=Append[{Abs@DegenBlock12[ToExpression@c,ToExpression@hl,ToExpression@hh,r,tL]}]@plotVector;
plotLegends=Append["Exact Degenerate \!\(\*SubscriptBox[\(h\), \(12\)]\)"]@plotLegends;
];
If[MemberQ[compareVec, "21"],
plotVector=Append[{Abs@DegenBlock21[ToExpression@c,ToExpression@hl,ToExpression@hh,r,tL]}]@plotVector;
plotLegends=Append["Exact Degenerate \!\(\*SubscriptBox[\(h\), \(21\)]\)"]@plotLegends;
];
Print[plotVector/.tL->1//N];
Print[Abs@SemiClassical[ToExpression@c,ToExpression@hl,ToExpression@hh,r,1]];
Print[LogPlot[Evaluate[plotVector],{tL,0.1,endTime},PlotRange->All,PlotLegends->plotLegends,PlotLabel->"Block in Lorentzian time",AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(L\)]\)","V(\!\(\*SubscriptBox[\(t\), \(L\)]\))"},PlotPoints->OptionValue[PointsPerTL]*endTime]];
,{i,startingRun,Length@results/2,runStep}];
];
End[]

EndPackage[]



