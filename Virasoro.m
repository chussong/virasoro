(* ::Package:: *)

BeginPackage["Virasoro`"]

(*For simple computations*)
Gethmn::usage = "Gethmn[m_,n_,b_] gives degenerate Virasoro operator dimension h_mn at the given b.";
Get\[Lambda]sq::usage = "Get\[Lambda]sq[p_,q_,b_] gives square of internal \[Lambda]_pq at the given b.";
Get\[Lambda]::usage = "Get\[Lambda][p_,q_,b_] gives internal \[Lambda]_pq at the given b.";
GetBsq::usage = "GetBsq[c_] gives b^2 corresponding to the given c.";
GetC::usage = "GetC[bsq_] gives c corresponding to the given b^2.";

(*For plotting in Lorentzian time*)
EKMono::usage = "EKMono[r_,tL_] gives the EllipticK[z] for these parameters accounting for the monodromy.";
qVal::usage = "qVal[r_,tL] gives the q value corresponding to the given Lorentzian time tL using radius r.";
VofT::usage = "QtoT[coeVec_,c_,hL_,hH_,hp_,r_,tL_,maxOrder_] gives the Lorentzian Virasoro block approximated by coeVec at Lorentzian time tL.";
SemiClassical::usage = "SemiClassical[c_,hL_,hH_,r_,tL_] gives the value of the semiclassical approximation of the vacuum block at Lorentzian time tL.";
DegenBlock12::usage = "DegenBlock12[c_,hL_,hh_,r_,tL_] gives the exact degenerate h_12 vacuum block at time tL.";
DegenBlock21::usage = "DegenBlock21[c_,hL_,hh_,r_,tL_] gives the exact degenerate h_21 vacuum block at time tL.";

(*For reading results from virasoro batch runs*)
VRun::usage = "VRun[c_, hl_, hh_, hp_, maxOrder_] or VRun[{c, hl, hh, hp, maxOrder}] computes the coefficients for the given parameter set(s) and returns them as a list of lists using the same format as VRead."
VRunFromFile::usage = "VRunFromFile[filename_] or simply VRun[filename_] invokes a C++ run from the named runfile and returns the results.";
VLengthCheck::usage = "VLengthCheck[filename_] gives the length of the result vector filename.";
VRead::usage = "VRead[filename_] gives the run results as a vector of vectors of numbers, ready to be used in Mathematica.";
VWrite::usage = "VWrite[results_,filename_] writes a results to a file with the given filename.";
VMakePlotLabel::usage = "VMakePlotLabel[results_,runNumber_] is an internal function that returns a list of parameters used in a run. This is put at the top of most plots.";
VPlotCoeffs::usage = "VPlotCoeffs[results_] gives a log-log plot of the raw coefficients. Options StartingRun, EndingRun, and RunStep specified as Option\[Rule]Value.";
VPlot::usage = "VPlot[results_(,Option\[Rule]Value)] plots Virasoro blocks from all runs in the VRead output results as a function of Lorentzian time, using radius r. runStep>1 allows you to skip runs.

Options are specified as in Wolfram's Plot[]: Option\[Rule]Value, e.g. EndTime\[Rule]40. Options are PlotScale (\"Linear\", \"SemiLog\", or \"LogLog\"), StartTime, EndTime, PointsPerTL, StartingRun, EndingRun, RunStep, and Compare. Setting Compare\[Rule]\"Semi\" compares to semiclassical vacuum block; other possibilities are \"12\" and \"21\" which compare to exact degenerate blocks.";
VConvByOrder::usage = "VConvByOrder[results_,tL_(,Option\[Rule]Value)] gives plots showing how adding more terms to each run in results_ improves the convergence at particular radius r_ and Lorentzian time tL_.

Options are specified as in Wolfram's Plot[]: Option\[Rule]Value, e.g. StartingRun\[Rule]10. Options are StartingRun, RunStep, and EndingRun.";
VConvByTL::usage = "VConvByTL[results_(,Option\[Rule]Value)] gives plots showing how the convergence of the last few orders of each run in results_ varies with tL at a given r_.

Options are specified as in Wolfram's Plot[]: Option\[Rule]Value, e.g. StartingRun\[Rule]10. Options are StartingRun, RunStep, EndingRun, StartTime, and EndTime.";

Begin["VirasoroInternal`"]
Gethmn[m_,n_,b_]:=b^2*(1-n^2)/4+(1-m^2)/(4*b^2)+(1-m*n)/2;
Get\[Lambda]sq[p_,q_,b_]:=4*Gethmn[p,q]-(b+1/b)^2;
Get\[Lambda][p_,q_,b_]:=Sqrt[Get\[Lambda]sq[p,q,b]];
GetBsq[c_]:=(-13+c+Sqrt[25-26c+c^2])/12;
GetC[bsq_]:=13+6*(bsq+1/bsq);

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
bsq=GetBsq[c];
1/(1-r*Exp[-I tL])^(2hL) r^b1 Exp[-I b1 tL/2]((Gamma[a1+b1-c1] Gamma[c1])/(Gamma[a1] Gamma[b1]) r^(-a1-b1+c1) Exp[-I (-tL a1-tL b1+tL c1)] Hypergeometric2F1[-a1+c1,-b1+c1,1-a1-b1+c1,r E^(-I tL)] +(Gamma[-a1-b1+c1] Gamma[c1])/(Gamma[c1-a1] Gamma[c1-b1]) Hypergeometric2F1[a1,b1,a1+b1-c1+1,r E^(-I tL)] )/.{a1->1+1/bsq,b1->(1+bsq+Sqrt[1+bsq^2+bsq *(2-4 hh)])/bsq,c1->2+2/bsq}/.z->1-r E^(-I tL)
];
DegenBlock21[c_,hL_,hh_,r_,tL_]:=Module[{a1,b1,c1,z,bsq},
bsq=GetBsq[c];
1/(1-r Exp[-I tL])^(2hL) r^b1 Exp[-I b1 tL/2]((Gamma[a1+b1-c1] Gamma[c1])/(Gamma[a1] Gamma[b1]) r^(-a1-b1+c1) Exp[-I (-tL a1-tL b1+tL c1)] Hypergeometric2F1[-a1+c1,-b1+c1,1-a1-b1+c1,r E^(-I tL)] +(Gamma[-a1-b1+c1] Gamma[c1])/(Gamma[c1-a1] Gamma[c1-b1]) Hypergeometric2F1[a1,b1,a1+b1-c1+1,r E^(-I tL)] )/.{a1->1+bsq,b1->(1+1/bsq+Sqrt[1+bsq^(-2)+1/bsq *(2-4 hh)])*bsq,c1->2+2*bsq}/.z->1-r E^(-I tL)
];

VRun::paramError = "Enter 5 parameters, either separately or as a vector.";
VRun::noVWSTP = "VWSTP not found. It must be in the same directory as Virasoro.m and must be marked as executable.";
VRun[c_,hl_,hh_,hp_,maxOrder_] := Module[{link,params,results},
link = Install[NotebookDirectory[]<>"vwstp"];
If[link==$Failed,Failure["C++NotFound",<|"MessageTemplate":>VRun::noVWSTP|>]];
params = ToString/@N[{c,hl,hh,hp,maxOrder},768];
Do[If[StringContainsQ["I"]@params[[i]],
params[[i]] = StringInsert[params[[i]],"(",1];
params[[i]] = StringInsert[params[[i]],")",-1];
params[[i]] = StringReplace[params[[i]], {" +" -> "", " -" -> "", " I" -> ""}];
],{i,1,Length@params}];
results = Global`VPass[params[[1]],params[[2]],params[[3]],params[[4]],params[[5]]];
Uninstall[link];
Return[results];
];
VRunFromFile::notFound = "The argument given was not a vector of parameters or valid filename.";
VRunFromFile[filename_]:=Module[{link,results},
If[FindFile[NotebookDirectory[]<>filename] == $Failed,
Failure["FileNotFound", <|"MessageTemplate":>VRunFromFile::notFound|>]];
link = Install[NotebookDirectory[]<>"vwstp"];
If[link==$Failed,Failure["C++NotFound",<|"MessageTemplate":>VRun::noVWSTP|>]];
results = Global`VPassFilename[NotebookDirectory[]<>filename];
Uninstall[link];
Return[results];
];
VRun[paramVec_]:=Module[{stringParam,results},
If[VectorQ[paramVec], stringParam=ToString/@paramVec, stringParam=StringReplace[ToString@paramVec," "->""]];
Switch[Length@stringParam
,0,results=VRunFromFile[stringParam];
,1,results=VRunFromFile[stringParam[[1]]];
,5,results=VRun[stringParam[[1]],stringParam[[2]],stringParam[[3]],stringParam[[4]],stringParam[[5]]];
,_, Failure["InvalidParameters", <|"MessageTemplate":>VRun::paramError|>]];
Return[results];
];
VLengthCheck[filename_]:=Module[{resultFile,entries, null},
resultFile=OpenRead[filename,BinaryFormat->True];
entries =0;
null[_String]:=Null;
entries=Length@ReadList[filename,null@String,NullRecords->True];
Close[filename];
Return[entries];
];
VRead[filename_]:=Module[{resultFile,entries,null,results},
entries=VLengthCheck[filename];
results=ConstantArray[0,entries];
resultFile=OpenRead[filename,BinaryFormat->True];
Do[
results[[2i-1]]=ToExpression@ReadLine[resultFile];
results[[2i]]=ToExpression@ReadLine[resultFile];
,{i,1,entries/2}];
Close[filename];
Return[results];
];
VWrite[results_, rawFilename_]:=Module[{filename, resultFile},
If[StringCount[rawFilename, "/"] == 0, filename = NotebookDirectory[]<>ToString@rawFilename, filename = ToString@rawFilename];
resultFile=OpenWrite[filename,PageWidth->Infinity];
Do[
WriteString[resultFile,NumberForm[results[[i]], ExponentFunction->(Null&)]];
WriteString[resultFile,"\n"];
,{i,1,Length@results}];
Close[resultFile];
Return[];
];
VMakePlotLabel[results_,runNumber_]:=Module[{c,hl,hh,hp,label},
c=results[[2runNumber-1]][[1]];
hl=results[[2runNumber-1]][[2]];
hh=results[[2runNumber-1]][[3]];
hp=results[[2runNumber-1]][[4]];
label="c="<>StringTake[ToString@c,Min[5,StringLength[ToString@c]]]<>"    \!\(\*SubscriptBox[\(h\), \(l\)]\)="<>StringTake[ToString@hl,Min[5,StringLength[ToString@hl]]]<>"    \!\(\*SubscriptBox[\(h\), \(h\)]\)="<>StringTake[ToString@hh,Min[5,StringLength[ToString@hh]]]<>"    \!\(\*SubscriptBox[\(h\), \(p\)]\)="<>StringTake[ToString@hp,Min[5,StringLength[ToString@hp]]];
Return[label];
];
Options[VPlotCoeffs]={StartingRun->1, EndingRun->0, RunStep->1};
VPlotCoeffs[results_,OptionsPattern[]]:=Module[{c,hl,hh,hp},
Do[If[Length@results[[2*i]]<10,Continue[]];
c=results[[2i-1]][[1]];
hl=results[[2i-1]][[2]];
hh=results[[2i-1]][[3]];
hp=results[[2i-1]][[4]];
Print[ListLogLogPlot[{results[[2*i]],-results[[2*i]]},PlotLabel->VMakePlotLabel[results,i],PlotMarkers->".",PlotStyle->{Lighter@Blue,Lighter@Red},PlotLegends->{"Blue > 0","Red < 0"},DataRange->{0,2*Length@results[[2*i]]}]];
(*Print[ListPlot[LogFluct[results[[2*i]]],PlotLabel\[Rule]"Fluctuations about smoothed average log", PlotMarkers\[Rule]".",ColorFunction\[Rule]Coloring,ColorFunctionScaling\[Rule]False,DataRange\[Rule]{0,2*Length@results[[2*i]]}]];*)
,{i,OptionValue[StartingRun],If[OptionValue[EndingRun]!=0,OptionValue[EndingRun],Length@results/2],OptionValue[RunStep]}];
];
Options[VPlot]={StartingRun->1, EndingRun->0, RunStep->1, r->0.99, PlotScale->"LogLog", StartTime->0.1, EndTime->30, Compare->{}, PointsPerTL->1};
VPlot[results_,OptionsPattern[]]:=Module[{c,hl,hh,hp,startTime,endTime, compareVec, plotVector, plotLegends},
startTime=OptionValue[StartTime];
endTime=OptionValue[EndTime];
If[VectorQ[OptionValue[Compare]],compareVec = OptionValue[Compare], compareVec = {OptionValue[Compare]}];
Do[If[Length@results[[2*i]]<10,Continue[]];
c=results[[2i-1]][[1]];
hl=results[[2i-1]][[2]];
hh=results[[2i-1]][[3]];
hp=results[[2i-1]][[4]];
plotVector := {Abs@VofT[results[[2*i]],c,hl,hh,hp,OptionValue[r],tL,2*Length@results[[2*i]]-2]};
plotLegends := {"Computed"};
If[MemberQ[compareVec, "Semi"],
plotVector=Append[Abs@SemiClassical[c,hl,hh,OptionValue[r],tL]]@plotVector;
plotLegends=Append["Semiclassical"]@plotLegends;
];
If[MemberQ[compareVec, "12"],
plotVector=Append[{Abs@DegenBlock12[c,hl,hh,OptionValue[r],tL]}]@plotVector;
plotLegends=Append["Exact Degenerate \!\(\*SubscriptBox[\(h\), \(12\)]\)"]@plotLegends;
];
If[MemberQ[compareVec, "21"],
plotVector=Append[{Abs@DegenBlock21[c,hl,hh,OptionValue[r],tL]}]@plotVector;
plotLegends=Append["Exact Degenerate \!\(\*SubscriptBox[\(h\), \(21\)]\)"]@plotLegends;
];
(*Print[plotVector/.tL\[Rule]1//N];
Print[Abs@SemiClassical[c,hl,hh,r,1]];*)
If[StringMatchQ[OptionValue[PlotScale],"Linear"],Print[Plot[Evaluate[plotVector],{tL,startTime,endTime},PlotRange->All,PlotLegends->plotLegends,PlotLabel->VMakePlotLabel[results,i],AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(L\)]\)","V(\!\(\*SubscriptBox[\(t\), \(L\)]\))"},PlotPoints->Max[OptionValue[PointsPerTL]*Ceiling[endTime-startTime],50]]]];
If[StringMatchQ[OptionValue[PlotScale],"SemiLog"]||StringMatchQ[OptionValue[PlotScale],"LogLinear"],Print[LogPlot[Evaluate[plotVector],{tL,startTime,endTime},PlotRange->All,PlotLegends->plotLegends,PlotLabel->VMakePlotLabel[results,i],AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(L\)]\)","V(\!\(\*SubscriptBox[\(t\), \(L\)]\))"},PlotPoints->Max[OptionValue[PointsPerTL]*Ceiling[endTime-startTime],50]]]];
If[StringMatchQ[OptionValue[PlotScale],"LogLog"],Print[LogLogPlot[Evaluate[plotVector],{tL,startTime,endTime},PlotRange->All,PlotLegends->plotLegends,PlotLabel->VMakePlotLabel[results,i],AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(L\)]\)","V(\!\(\*SubscriptBox[\(t\), \(L\)]\))"},PlotPoints->Max[OptionValue[PointsPerTL]*Ceiling[endTime-startTime],50]]]];
,{i,OptionValue[StartingRun],If[OptionValue[EndingRun]!=0,OptionValue[EndingRun],Length@results/2],OptionValue[RunStep]}];
];
Options[VConvByOrder]={r->0.99, StartingRun->1, RunStep->1, EndingRun->0};
VConvByOrder[results_,tL_,OptionsPattern[]]:=Module[{plotLabel,q},
Do[
Print[DiscretePlot[Table[q^k,{k,0,2*Floor[i/2],2}].Take[results[[entry]],Floor[i/2]+1]/.q->qVal[OptionValue[r],tL]//Abs,{i,Length[results[[entry]]]/5,2*Length[results[[entry]]],2},PlotRange->Full,AxesLabel->{"Max Order","H(r="<>ToString@OptionValue[r]<>",tL="<>ToString@tL<>")"},Filling->None,PlotLabel->VMakePlotLabel[results,entry/2]]];
,{entry,2*OptionValue[StartingRun],If[OptionValue[EndingRun]!=0,2*OptionValue[EndingRun],Length@results],2*OptionValue[RunStep]}];
];
Options[VConvByTL]={StartingRun->1, RunStep->1, EndingRun->0,StartTime->0.1,EndTime->50,r->0.99};
VConvByTL[results_,OptionsPattern[]]:=Module[{plotLabel,q},
Do[
Print[Plot[(Table[q^k,{k,0,2*Length@results[[entry]]-2,2}].results[[entry]])/(Table[q^k,{k,0,2*Max[Length@results[[entry]]-10,1]-2,2}].Take[results[[entry]],Max[Length@results[[entry]]-10,1]])/.q->qVal[OptionValue[r],tL]//Abs,{tL,OptionValue[StartTime],OptionValue[EndTime]},PlotRange->Full,AxesLabel->{"tL","H(r="<>ToString@OptionValue[r]<>")"},PlotLabel->VMakePlotLabel[results,entry/2]]];
,{entry,2*OptionValue[StartingRun],If[OptionValue[EndingRun]!=0,2*OptionValue[EndingRun],Length@results],2*OptionValue[RunStep]}];
];
End[]

EndPackage[]



