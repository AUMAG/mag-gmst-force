(* ::Package:: *)

(* ::Section:: *)
(*Citation for this work *)


(* ::Text:: *)
(*M. Forbes, W.S.P Robertson, A.C. Zander, J.J.H. Paulides (2024) "The Magnetic Field from Cylindrical Arc Coils and Magnets: A Compendium with New Analytic Solutions for Radial Magnetisation and Azimuthal Current" Advanced Physics Research vol. 3, no. 7, p. 2300136.*)


(* ::Text:: *)
(*@Article {Forbes2024,*)
(*        author = {Forbes, M. and Robertson, W.S.P. and Zander, A.C. and Paulides, J.J.H},*)
(*        journal = {Advanced Physics Research},*)
(*        title = {The Magnetic Field from Cylindrical Arc Coils and Magnets: A Compendium with New Analytic Solutions for Radial Magnetisation and Azimuthal Current},*)
(*        doi = {10.1002/apxr.202300136},*)
(*        publisher = {Wiley},*)
(*        number = {7},*)
(*        pages = {2300136},*)
(*        volume = {3},*)
(* }*)


(* ::Section:: *)
(**)


BeginPackage["MagCylField`"]


CylB\[ScriptD]::usage="Diametric magnetisation\n CylB\[ScriptD][M,\[CurlyPhi]\:2606,{\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[Rho]\), \(2\)]\)'},\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],{\!\(\*SubscriptBox[\(z\), \(1\)]\)',\!\(\*SubscriptBox[\(z\), \(2\)]\)'},z] with: M magnetisation magnitude (A/m), \[CurlyPhi]\:2606 magnetision direction (rad), r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the magnet surface are (potentially) undefined."
CylB\[Rho]::usage="Radial magnetisation\n CylB\[Rho][M,P,{\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[Rho]\), \(2\)]\)'},\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],{\!\(\*SubscriptBox[\(z\), \(1\)]\)',\!\(\*SubscriptBox[\(z\), \(2\)]\)'},z] with: M magnetisation magnitude (A/m), P number of series terms, r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the magnet surface are (potentially) undefined."
CylB\[CurlyPhi]::usage="Diametric magnetisation\n CylB\[CurlyPhi][M,{\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[Rho]\), \(2\)]\)'},\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],{\!\(\*SubscriptBox[\(z\), \(1\)]\)',\!\(\*SubscriptBox[\(z\), \(2\)]\)'},z] with: M magnetisation magnitude (A/m), r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the magnet surface are (potentially) undefined."
CylB\[ScriptZ]::usage="Diametric magnetisation\n CylB\[ScriptZ][M,{\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[Rho]\), \(2\)]\)'},\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],{\!\(\*SubscriptBox[\(z\), \(1\)]\)',\!\(\*SubscriptBox[\(z\), \(2\)]\)'},z] with: M magnetisation magnitude (A/m), r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the magnet surface are (potentially) undefined."
CylB\[ScriptCapitalI]\[ScriptF]::usage="Azmimuthal current filament\n CylB\[ScriptCapitalI]\[ScriptF][I,\[Rho]',\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],z',z] with: I current (A), r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the filament are undefined."
CylB\[ScriptCapitalK]\[ScriptD]::usage="Azmimuthal current disc\n CylB\[ScriptCapitalK]\[ScriptD][K,P,{\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[Rho]\), \(2\)]\)'},\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],z',z] with: K line current density (A/m), P number of series terms, r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the disc surface are undefined."
CylB\[ScriptCapitalK]\[ScriptS]::usage="Azmimuthal current shell\n CylB\[ScriptCapitalK]\[ScriptS][K,\[Rho]',\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],{\!\(\*SubscriptBox[\(z\), \(1\)]\)',\!\(\*SubscriptBox[\(z\), \(2\)]\)'},z] with: K line current density (A/m), r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the shell surface are undefined."
CylB\[ScriptCapitalJ]\[ScriptV]::usage="Azmimuthal current volume\n CylB\[ScriptCapitalJ]\[ScriptV][J,P,{\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[Rho]\), \(2\)]\)'},\[Rho],{\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(1\)]\)',\!\(\*SubscriptBox[\(\[CurlyPhi]\), \(2\)]\)'},\[CurlyPhi],{\!\(\*SubscriptBox[\(z\), \(1\)]\)',\!\(\*SubscriptBox[\(z\), \(2\)]\)'},z] with: J surface current density (A/\!\(\*SuperscriptBox[\(m\), \(2\)]\)), P number of series terms, r' source limits (m,rad,m), r field point (m,rad,m).\n Returns the magnetic field in cylindrical coordinates for \[Rho]\[NotEqual]0 and Cartesian coordinates for \[Rho]=0. \n Points on the coil surface are undefined. \n This is a computationally expensive calculation for high P."


Begin["`Private`"]
(*Needs["Carlson`"]*)


(* ::Section::Closed:: *)
(*Nomenclature*)


(*Constants*)
u0 = 4\[Pi]*10^-7
(*Commonly occuring variables/functions*)
\[CurlyRho][\[Rho]_,\[Rho]p_] := (\[Rho]+\[Rho]p)
\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho]_,\[Rho]p_] := (\[Rho]-\[Rho]p)
\[CapitalPhi][\[CurlyPhi]_,\[CurlyPhi]p_] := (\[CurlyPhi]-\[CurlyPhi]p)
Z[z_,zp_] := (z-zp)
\[Phi][\[CurlyPhi]_,\[CurlyPhi]p_] := 1/2 (\[CurlyPhi]-\[CurlyPhi]p+\[Pi])
k[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[((4\[Rho] \[Rho]p)/(R[\[Rho],\[Rho]p,z,zp]^2)) ]
\!\(\*OverscriptBox[\(k\), \(_\)]\)[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[(4\[Rho] \[Rho]p)/\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2] 
a[\[Rho]_,z_,zp_] := Sqrt[(2\[Rho])/(\[Rho]+L[\[Rho],z,zp])]
\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho]_,z_,zp_] := Sqrt[(2\[Rho])/(\[Rho]-L[\[Rho],z,zp])]
\[Kappa][\[Rho]_,\[Rho]p_]:= Sqrt[(4\[Rho] \[Rho]p)/\[CurlyRho][\[Rho],\[Rho]p]^2]
L[\[Rho]_,z_,zp_] := Sqrt[Z[z,zp]^2+\[Rho]^2]
R[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[Z[z,zp]^2+\[CurlyRho][\[Rho],\[Rho]p]^2]
\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[Z[z,zp]^2+\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]^2]
S[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[L[\[Rho],z,zp]+\[Rho]p]
\!\(\*OverscriptBox[\(S\), \(_\)]\)[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[L[\[Rho],z,zp]-\[Rho]p]
T[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[L[\[Rho],z,zp]^2+\[Rho]p^2]
\!\(\*OverscriptBox[\(T\), \(_\)]\)[\[Rho]_,\[Rho]p_,z_,zp_] := Sqrt[L[\[Rho],z,zp]^2-\[Rho]p^2]
G[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_]:=1/Sqrt[T[\[Rho],\[Rho]p,z,zp]^2-2\[Rho] \[Rho]p Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]]
ns[\[CapitalPhi]p_] := Floor[\!\(\*OverscriptBox[\(\[CapitalPhi]\), \(~\)]\)[\[CapitalPhi]p]]
qs[\[CapitalPhi]p_] := (-1)^Floor[2\!\(\*OverscriptBox[\(\[CapitalPhi]\), \(~\)]\)[\[CapitalPhi]p]+1]
ts[\[CapitalPhi]p_] := (1-2Abs[Round[\!\(\*OverscriptBox[\(\[CapitalPhi]\), \(~\)]\)[\[CapitalPhi]p]]-\!\(\*OverscriptBox[\(\[CapitalPhi]\), \(~\)]\)[\[CapitalPhi]p]])
zs[\[Rho]_,\[Rho]p_,z_,zp_]:= Sqrt[\[Rho]p^2/T[\[Rho],\[Rho]p,z,zp]^2]
ys[\[Rho]_,\[Rho]p_,z_,zp_]:= Sqrt[Z[z,zp]^2/T[\[Rho],\[Rho]p,z,zp]^2]
\[Omega][\[CurlyPhi]_,\[CurlyPhi]p_]:=(1+Round[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]/\[Pi]-1])
v[\[Nu]_] := (\[Nu]-2Floor[\[Nu]/2])
\!\(\*OverscriptBox[\(\[CapitalPhi]\), \(~\)]\)[\[CapitalPhi]p_] := \[CapitalPhi]p/(2\[Pi])
\[CapitalUpsilon][\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_] := (Z[z,zp] (\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]-\[Rho]p) )/(\[Rho] Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1)
(*transformed Legendre functions*)
EllipticFT[\[Phi]_,k_]:= 2ns[2\[Phi]-\[Pi]] EllipticK[k]+qs[2\[Phi]-\[Pi]]EllipticF[\[Pi]/2 ts[2\[Phi]-\[Pi]],k]
EllipticET[\[Phi]_,k_]:= 2ns[2\[Phi]-\[Pi]]EllipticE[k]+qs[2\[Phi]-\[Pi]]EllipticE[\[Pi]/2 ts[2\[Phi]-\[Pi]],k]
EllipticPiT[\[Alpha]_,\[Phi]_,k_]:= 2ns[2\[Phi]-\[Pi]]EllipticPi[\[Alpha],k]+qs[2\[Phi]-\[Pi]]EllipticPi[\[Alpha],\[Pi]/2 ts[2\[Phi]-\[Pi]],k]
EllipticDT[\[Phi]_,k_]:= 1/k (EllipticFT[\[Phi],k]-EllipticET[\[Phi],k])
EllipticD[k_]:= 1/k (EllipticK[k]-EllipticE[k])
(*Canonical regularised beta summands*)
\[CapitalXi][\[Rho]_,\[Rho]p_,z_,zp_,\[Nu]_,p_]:= 1/(1+2 p+\[Nu]) (\[Rho]/(2 L[\[Rho],z,zp]))^(1+2 p+\[Nu]) BetaRegularized[zs[\[Rho],\[Rho]p,z,zp]^2,1+p+\[Nu]/2,1/2+p+\[Nu]/2] Binomial[1+2 p+\[Nu],p]
\[Xi][\[Rho]_,\[Rho]p_,z_,zp_,\[Nu]_,p_]:=((\[Rho] \[Rho]p)/(\[Rho]^2+\[Rho]p^2))^(2 p+\[Nu]+1) BetaRegularized[ys[\[Rho],\[Rho]p,z,zp]^2,1/2,2p+\[Nu]+1] 1/(\[Nu]+2p+1) Binomial[2p+\[Nu]+1,p]
\[Tau][\[Rho]_,\[Rho]p_,z_,zp_,\[Nu]_,p_]:= Pochhammer[5/2+p+\[Nu],p+1/2]  (\[Rho]p/\[Rho])^\[Nu] \[Xi][\[Rho]p,\[Rho],z,zp,\[Nu]+1+2p,0] 
(*Elementary functions*)
\[Alpha]1[\[Rho]_,\[Rho]p_,z_,zp_] := ArcTanh[zs[\[Rho],\[Rho]p,z,zp]]
\[Alpha]2[\[Rho]_,\[Rho]p_,z_,zp_] := ArcTanh[ys[\[Rho],\[Rho]p,z,zp]]
\[Alpha]3[\[Rho]_,\[Rho]p_,z_,zp_]:= \[Rho] ArcTan[(Z[z,zp] \[Rho]p)/(\[Rho] T[\[Rho],\[Rho]p,z,zp])]-Z[z,zp] ArcTanh[\[Rho]p/T[\[Rho],\[Rho]p,z,zp]] - \[Rho]p Log[Z[z,zp]+T[\[Rho],\[Rho]p,z,zp]]
(*Single series*)
\[Beta]1[\[Rho]_,\[Rho]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\ \(\[CapitalXi][\[Rho], \[Rho]p, z, zp, \(-1\), p + 1]\)\)
\[Beta]2[\[Rho]_,\[Rho]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\ \(\[Xi][\[Rho], \[Rho]p, z, zp, \(-1\), p + 1]\)\)
\[Beta]3[\[Rho]_,\[Rho]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\(
\*FractionBox[\(\(\ \)\(Binomial[2 \((1 + p)\), 1 + p]\)\(\ \)\), \(
\*SuperscriptBox[\(4\), \(p\)] \(\((1 + p)\)!\)\)] \(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\[Tau][\[Rho], \[Rho]p, z, zp, \[Nu], p]\)\)\)
(*Double series*)
\[Gamma]1[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_] := 2\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\(Cos[\((\[Nu] + 1)\) \[CapitalPhi][\[CurlyPhi], \[CurlyPhi]p]] \(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\[CapitalXi][\[Rho], \[Rho]p, z, zp, \[Nu], p]\)\)\)
\[Gamma]2[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_] :=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\(Cos[\((\[Nu] + 1)\) \[CapitalPhi][\[CurlyPhi], \[CurlyPhi]p]] \(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\[Xi][\[Rho], \[Rho]p, z, zp, \[Nu], p]\)\)\)
\[Gamma]3[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_] := \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\(
\*FractionBox[\(\[Chi][\[CurlyPhi], \[CurlyPhi]p, p]\), \(Pochhammer[2 + 
\*FractionBox[\(p\), \(2\)], 1 + 
\*FractionBox[\(p\), \(2\)]]\)] \(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\[Tau][\[Rho], \[Rho]p, z, zp, \[Nu], 
\*FractionBox[\(p - 1\), \(2\)]]\)\)\)
(*Reduced form & computational series*)
\[Delta]1[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\ \(\[CapitalXi][\[Rho], \[Rho]p, z, zp, \[Nu], 0] \[Chi][\[CurlyPhi], \[CurlyPhi]p, \[Nu]]\)\)
\[Delta]2[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\ \(\[Xi][\[Rho], \[Rho]p, z, zp, \[Nu], 0] \[Chi][\[CurlyPhi], \[CurlyPhi]p, \[Nu]]\)\)
\[Delta]3[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\ \(\[CapitalXi][\[Rho], \[Rho]p, z, zp, \[Nu], 0] Sum[\[Lambda][\[CurlyPhi], \[CurlyPhi]p, \[Nu], p]\ Binomial[1 + \[Nu], p], {p, 0, Floor[\[Nu]/2]}]\)\)
\[Chi][\[CurlyPhi]_,\[CurlyPhi]p_,\[Nu]_]:=(2^\[Nu] (v[\[Nu]+1]Sign[Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]]+ v[\[Nu]]s[\[CurlyPhi],\[CurlyPhi]p]) Beta[Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]^2,1/2,1+\[Nu]/2]+v[\[Nu]](2^(\[Nu]+1) Beta[1/2,1+\[Nu]/2]\[Omega][\[CurlyPhi],\[CurlyPhi]p]-Binomial[1+\[Nu],(1+\[Nu])/2]\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]))
\[Chi]\[Zeta][\[CurlyPhi]_,\[CurlyPhi]p_,\[Nu]_] := 2^\[Nu] Sign[Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]]Beta[Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]^2,1/2,\[Nu]/2+1]
\[Chi]\[Eta][\[CurlyPhi]_,\[CurlyPhi]p_,\[Nu]_] := 2^\[Nu] s[\[CurlyPhi],\[CurlyPhi]p]Beta[Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]^2,1/2,\[Nu]/2+1]
\[Chi]\[Iota][\[CurlyPhi]_,\[CurlyPhi]p_,\[Nu]_] := 2^(\[Nu]+1) Beta[1/2,\[Nu]/2+1]\[Omega][\[CurlyPhi],\[CurlyPhi]p]-Binomial[\[Nu]+1,\[Nu]/2+1/2]\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]
\[Zeta]1[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,p_] :=  \[CapitalXi][\[Rho],\[Rho]p,z,zp,2p,0]\[Chi]\[Zeta][\[CurlyPhi],\[CurlyPhi]p,2p]
\[Eta]1[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,p_] :=  \[CapitalXi][\[Rho],\[Rho]p,z,zp,2p+1,0]\[Chi]\[Eta][\[CurlyPhi],\[CurlyPhi]p,2p+1]
\[Iota]1[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,p_] :=  \[CapitalXi][\[Rho],\[Rho]p,z,zp,2p+1,0]\[Chi]\[Iota][\[CurlyPhi],\[CurlyPhi]p,2p+1]
\[Zeta]2[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,p_] :=  \[Xi][\[Rho],\[Rho]p,z,zp,2p,0]\[Chi]\[Zeta][\[CurlyPhi],\[CurlyPhi]p,2p]
\[Eta]2[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,p_] :=  \[Xi][\[Rho],\[Rho]p,z,zp,2p+1,0]\[Chi]\[Eta][\[CurlyPhi],\[CurlyPhi]p,2p+1]
\[Iota]2[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,p_] :=  \[Xi][\[Rho],\[Rho]p,z,zp,2p+1,0]\[Chi]\[Iota][\[CurlyPhi],\[CurlyPhi]p,2p+1]
\[Zeta]3[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,\[Nu]_,p_] := \[Tau][\[Rho],\[Rho]p,z,zp,\[Nu],p-1/2] \[Chi]\[Zeta][\[CurlyPhi],\[CurlyPhi]p,2p]/Pochhammer[2+p,1+p]
\[Eta]3[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,\[Nu]_,p_] := \[Tau][\[Rho],\[Rho]p,z,zp,\[Nu],p] \[Chi]\[Eta][\[CurlyPhi],\[CurlyPhi]p,2p+1]/Pochhammer[5/2+p,3/2+p]
\[Iota]3[\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,\[Nu]_,p_] := \[Tau][\[Rho],\[Rho]p,z,zp,\[Nu],p] \[Chi]\[Iota][\[CurlyPhi],\[CurlyPhi]p,2p+1]/Pochhammer[5/2+p,3/2+p]
\[Delta]1\[Zeta]\[Eta]\[Iota][\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\((\[Zeta]1[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, p] + \[Eta]1[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, p] + \[Iota]1[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, p])\)\)
\[Delta]2\[Zeta]\[Eta]\[Iota][\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\((\[Zeta]2[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, p] + \[Eta]2[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, p] + \[Iota]2[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, p])\)\)
\[Gamma]3\[Zeta]\[Eta]\[Iota][\[Rho]_,\[Rho]p_,\[CurlyPhi]_,\[CurlyPhi]p_,z_,zp_,P_]:= \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(P\)]\(
\*UnderoverscriptBox[\(\[Sum]\), \(p = 0\), \(P\)]\((\[Zeta]3[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, \[Nu], p] + \[Eta]3[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, \[Nu], p] + \[Iota]3[\[Rho], \[Rho]p, \[CurlyPhi], \[CurlyPhi]p, z, zp, \[Nu], p])\)\)\)
(*Misc*)
s[\[CurlyPhi]_,\[CurlyPhi]p_]:=If[Sin[2 \[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]==0,(-1)^Round[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]/\[Pi]+1/2], Sign[Sin[2\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]]]
\[Lambda][\[CurlyPhi]_,\[CurlyPhi]p_,\[Nu]_,p_]:=If[\[Nu]-2p==1,\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p],Sin[(-1+\[Nu]-2p)\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]/(-1+\[Nu]-2p)]+(2 Sin[(1+\[Nu]-2p) \[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]])/(1+\[Nu]-2p)+Sin[(3+\[Nu]-2p)\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]/(3+\[Nu]-2p)
\[CurlyTheta][\[Rho]_,\[Rho]p_,z_,zp_,P_] :=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(\[Nu] = 0\), \(Floor[
\*FractionBox[\(P - 1\), \(2\)]]\)]\ \(\[CapitalXi][\[Rho]p, \[Rho], z, zp, 2  \[Nu] + 1, 0]\ Binomial[2 + 2  \[Nu], \[Nu]]\)\)
(*check if field point within volume*)
InsideVolume[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] :=  (\[Rho] <\[Rho]p[[2]]) && (\[Rho] > \[Rho]p[[1]]) && (z <zp[[2]]) && (z> zp[[1]]) && If[Mod[\[CurlyPhi]p[[2]]-\[CurlyPhi]p[[1]],2\[Pi]] < \[Pi],Mod[\[CurlyPhi]-\[CurlyPhi]p[[1]],2\[Pi]]<Mod[\[CurlyPhi]p[[2]]-\[CurlyPhi]p[[1]],2\[Pi]],Mod[\[CurlyPhi]p[[2]]-\[CurlyPhi]p[[1]],2\[Pi]]<Mod[\[CurlyPhi]-\[CurlyPhi]p[[1]],2\[Pi]]]
InsideVolumeAxis[\[Rho]p_,zp_,z_] := (0 == \[Rho]p[[1]]) && (z < zp[[2]]) && (z > zp[[1]]) 


(* ::Section::Closed:: *)
(*Diametric magnetisation*)


(*Rename and format result*)
CylB\[ScriptD][M_,\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[BdAna[M,\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
BdAna[M_,\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BdAnaAxis[M,\[CurlyPhi]\:2606,\[Rho]p,\[CurlyPhi]p,zp,z] , (M u0)/(4\[Pi]) (Sum[(-1)^(m+n+q) BdSummand[\[CurlyPhi]\:2606,\[Rho]p[[m]],\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z],{m,1,2},{n,1,2},{q,1,2}]+ Md[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z])]
BdAnaAxis[M_,\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (M u0)/(4\[Pi]) (Sum[(-1)^(m+n+q) BdSummandAxis[\[CurlyPhi]\:2606,\[Rho]p[[m]],\[CurlyPhi]p[[q]],zp[[n]],z],{m,1,2},{n,1,2},{q,1,2}] + MdAxis[\[CurlyPhi]\:2606,\[Rho]p,zp,z])

(*Magnetisation vector*)
Md[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[InsideVolume[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],4\[Pi]{Cos[\[CurlyPhi]\:2606-\[CurlyPhi]],Sin[\[CurlyPhi]\:2606-\[CurlyPhi]],0},{0,0,0}]

MdAxis[\[CurlyPhi]\:2606_,\[Rho]p_,zp_,z_] := If[InsideVolumeAxis[\[Rho]p,zp,z],4\[Pi]{Cos[\[CurlyPhi]\:2606],Sin[\[CurlyPhi]\:2606],0},{0,0,0}]

(*Geometry special cases*)
BdAna[M_,\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BdAnaAxis[M,\[CurlyPhi]\:2606,\[Rho]p,{0,2\[Pi]},zp,z],(M u0)/(4\[Pi]) (Sum[(-1)^(m+n ) BdSummandAS[\[CurlyPhi]\:2606,\[Rho]p[[m]],\[Rho],\[CurlyPhi],zp[[n]],z],{m,1,2},{n,1,2}] + Md[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],{0,2\[Pi]},\[Pi],zp,z])]
BdAnaAxis[M_,\[CurlyPhi]\:2606_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (M u0)/(4\[Pi]) (Sum[(-1)^(m+n) BdSummandAxis[\[CurlyPhi]\:2606,\[Rho]p[[m]],{0,2\[Pi]},zp[[n]],z],{m,1,2},{n,1,2}] + MdAxis[\[CurlyPhi]\:2606,\[Rho]p,zp,z])
BdAna[M_,\[CurlyPhi]\:2606_,{0,\[Rho]p_},\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := (M u0)/(4\[Pi]) (Sum[(-1)^(n+q) (BdSummand1[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z] + Sum[(-1)^m  BdSummand2[\[CurlyPhi]\:2606,(m-1)\[Rho]p,\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z],{m,1,2}]),{n,1,2},{q,1,2}]+ Md[\[CurlyPhi]\:2606,{0,\[Rho]p},\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z])
BdAna[M_,\[CurlyPhi]\:2606_,{0,\[Rho]p_},\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := (M u0)/(4\[Pi]) (Sum[(-1)^n  BdSummandAS[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi],zp[[n]],z],{n,1,2}] + Md[\[CurlyPhi]\:2606,{0,\[Rho]p},\[Rho],{0,2\[Pi]},\[Pi],zp,z])
(*Summands*)
BdSummand[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BdS\[Rho]1[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + BdS\[Rho]2[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BdS\[CurlyPhi]1[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + BdS\[CurlyPhi]2[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BdSz1[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + BdSz2[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}
BdSummand1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BdS\[Rho]1[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BdS\[CurlyPhi]1[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BdSz1[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}
BdSummand2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BdS\[Rho]2[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BdS\[CurlyPhi]2[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BdSz2[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}

BdS\[Rho]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := 1/\[Rho] (\[Rho]p Z[z,zp])/ R[\[Rho],\[Rho]p,z,zp] Cos[\[CurlyPhi]\:2606-\[CurlyPhi]](2EllipticDT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-(1+(2 \!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p])/(\[Kappa][\[Rho], \[Rho]p]^2 \[CurlyRho][\[Rho],\[Rho]p]))EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-(\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p] (\[Kappa][\[Rho],\[Rho]p]^2-2))/(\[CurlyRho][\[Rho],\[Rho]p] \[Kappa][\[Rho],\[Rho]p]^2) EllipticPiT[\[Kappa][\[Rho],\[Rho]p]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]) +1 /(2\[Rho]^2) Sin[\[CurlyPhi]\:2606-\[CurlyPhi]](-Z[z,zp] G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1+\[CurlyRho][\[Rho],\[Rho]p]\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p] ArcTanh[G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1/Z[z,zp]]) 
BdS\[Rho]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Which[(\[CurlyPhi]p==\[CurlyPhi]+\[Pi]) || (\[CurlyPhi]p==\[CurlyPhi]-\[Pi]),ArcTanh[R[\[Rho],\[Rho]p,z,zp]/Z[z,zp]] Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p],True,Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p](ArcTan[\[CapitalUpsilon][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]]Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]-Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]] ArcTanh[G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1/Z[z,zp]])]
BdS\[CurlyPhi]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := 2/\[Rho] \[Rho]p Z[z,zp]  Sin[\[CurlyPhi]\:2606-\[CurlyPhi]]/R[\[Rho],\[Rho]p,z,zp] (\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]^2/(4 \[Rho] \[Rho]p) (EllipticPiT[\[Kappa][\[Rho],\[Rho]p]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2])-EllipticDT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2])-1/(2 \[Rho]^2) Cos[\[CurlyPhi]\:2606-\[CurlyPhi]](Z[z,zp]G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1+(\[Rho]^2+\[Rho]p^2) ArcTanh[G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1/Z[z,zp]])  
BdS\[CurlyPhi]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Which[(\[CurlyPhi]p==\[CurlyPhi]+\[Pi])|| (\[CurlyPhi]p==\[CurlyPhi]-\[Pi]),0,True,Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p](ArcTan[\[CapitalUpsilon][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]]Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]] + Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]ArcTanh[G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1/Z[z,zp]])]
BdSz1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := (2 \[Rho]p Cos[\[CurlyPhi]\:2606-\[CurlyPhi]])/R[\[Rho],\[Rho]p,z,zp] (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-2EllipticDT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2])+ Sin[\[CurlyPhi]\:2606-\[CurlyPhi]]/ \[Rho] G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1 
BdSz2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Log[\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]+G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1] Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]

BdSummandAxis[\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {BdSAxisx[\[CurlyPhi]\:2606,\[Rho]p,\[CurlyPhi]p,zp,z], BdSAxisy[\[CurlyPhi]\:2606,\[Rho]p,\[CurlyPhi]p,zp,z], BdSAxisz[\[CurlyPhi]\:2606,\[Rho]p,\[CurlyPhi]p,zp,z]}

BdSAxisx[\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] :=  (Z[z,zp] (2\[CurlyPhi]p Cos[\[CurlyPhi]\:2606]-Sin[\[CurlyPhi]\:2606-2\[CurlyPhi]p]))/(4 Sqrt[Z[z,zp]^2+\[Rho]p^2])-ArcTanh[Sqrt[Z[z,zp]^2+\[Rho]p^2]/Z[z,zp]]Cos[\[CurlyPhi]p]Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]
BdSAxisy[\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (Z[z,zp] (-Cos[\[CurlyPhi]\:2606-2 \[CurlyPhi]p]+2 \[CurlyPhi]p Sin[\[CurlyPhi]\:2606]))/(4 Sqrt[Z[z,zp]^2+\[Rho]p^2])-ArcTanh[Sqrt[Z[z,zp]^2+\[Rho]p^2]/Z[z,zp]]Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]Sin[\[CurlyPhi]p]
BdSAxisz[\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (-(\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2])+ArcTanh[\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2]]) Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]


BdSummandAS[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := {BdSAS\[Rho][\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BdSAS\[CurlyPhi][\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BdSASz[\[CurlyPhi]\:2606,\[Rho]p,\[Rho],\[CurlyPhi],zp,z]}

BdSAS\[Rho][\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 2Cos[\[CurlyPhi]\:2606-\[CurlyPhi]] (\[Rho]p Z[z,zp])/( \[Rho] R[\[Rho],\[Rho]p,z,zp]) ((1+(2 \!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p])/(\[Kappa][\[Rho], \[Rho]p]^2 \[CurlyRho][\[Rho],\[Rho]p]))EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]-2EllipticD[k[\[Rho],\[Rho]p,z,zp]^2]+(\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p] (\[Kappa][\[Rho],\[Rho]p]^2-2))/(\[CurlyRho][\[Rho],\[Rho]p] (\[Kappa][\[Rho],\[Rho]p]^2) ) EllipticPi[\[Kappa][\[Rho],\[Rho]p]^2,k[\[Rho],\[Rho]p,z,zp]^2]) 
BdSAS\[CurlyPhi][\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 4 Sin[\[CurlyPhi]\:2606-\[CurlyPhi]] (   \[Rho]p Z[z,zp])/( \[Rho] R[\[Rho],\[Rho]p,z,zp]) (EllipticD[k[\[Rho],\[Rho]p,z,zp]^2]+\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]^2/(4 \[Rho] \[Rho]p) (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]-EllipticPi[\[Kappa][\[Rho],\[Rho]p]^2,k[\[Rho],\[Rho]p,z,zp]^2]))
BdSASz[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := (4 \[Rho]p Cos[\[CurlyPhi]\:2606-\[CurlyPhi]])/R[\[Rho],\[Rho]p,z,zp] (2EllipticD[k[\[Rho],\[Rho]p,z,zp]^2]-EllipticK[k[\[Rho],\[Rho]p,z,zp]^2])
(*Along the axis*)
BdSAxisx[\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,zp_]:=0
BdSAxisy[\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,zp_]:=0
BdSAxisz[\[CurlyPhi]\:2606_,\[Rho]p_,\[CurlyPhi]p_,zp_,zp_]:=(Log[\[Rho]p]-1) Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]
(*Along the axis & axisymmetric*)
BdSAxisx[\[CurlyPhi]\:2606_,\[Rho]p_,{0,2\[Pi]},zp_,z_] :=  (\[Pi] Z[z,zp] Cos[\[CurlyPhi]\:2606])/Sqrt[Z[z,zp]^2+\[Rho]p^2]
BdSAxisy[\[CurlyPhi]\:2606_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (\[Pi] Z[z,zp] Sin[\[CurlyPhi]\:2606])/Sqrt[Z[z,zp]^2+\[Rho]p^2]
BdSAxisz[\[CurlyPhi]\:2606_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BdSAxisx[\[CurlyPhi]\:2606_,\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
BdSAxisy[\[CurlyPhi]\:2606_,\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
BdSAxisz[\[CurlyPhi]\:2606_,\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
(*On the shell plane*)
EllipticPiT[1,\[Phi]_,k_] := EllipticFT[\[Phi],k]-1/(1-k)(EllipticET[\[Phi],k]-Sqrt[1-k Sin[\[Phi]]^2] Tan[\[Phi]]);
(*On the shell plane & axisymmetric*)
BdSAS\[Rho][\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]_,zp_,z_] := 2 ( Z[z,zp] Cos[\[CurlyPhi]\:2606-\[CurlyPhi]] )/ R[\[Rho]p,\[Rho]p,z,zp]  (EllipticK[k[\[Rho]p,\[Rho]p,z,zp]^2]-2 EllipticD[k[\[Rho]p,\[Rho]p,z,zp]^2])
BdSAS\[CurlyPhi][\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]_,zp_,z_] := 4 Z[z,zp]  Sin[\[CurlyPhi]\:2606-\[CurlyPhi]]/R[\[Rho]p,\[Rho]p,z,zp] EllipticD[k[\[Rho]p,\[Rho]p,z,zp]^2]
(*On the section plane*)
BdS\[Rho]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := -ArcTanh[\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]/Z[z,zp]] Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]
BdS\[CurlyPhi]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
(*On the disc plane*)
BdS\[Rho]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BdS\[Rho]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
BdS\[CurlyPhi]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BdS\[CurlyPhi]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the axial line*)
BdS\[Rho]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := Cos[\[CurlyPhi]\:2606-\[CurlyPhi]p]  Z[z,zp]/ R[\[Rho]p,\[Rho]p,z,zp] (EllipticK[k[\[Rho]p,\[Rho]p,z,zp]^2]-2EllipticD[k[\[Rho]p,\[Rho]p,z,zp]^2])-1/(2 \[Rho]p^2) Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]Sign[Z[z,zp]]Z[z,zp]^2
BdS\[Rho]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := -Log[Abs[Z[z,zp]]]Sign[Z[z,zp]]Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]
BdS\[CurlyPhi]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 2 Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]    Z[z,zp]/ R[\[Rho]p,\[Rho]p,z,zp] EllipticD[k[\[Rho]p,\[Rho]p,z,zp]^2]-Sign[Z[z,zp]]Cos[\[CurlyPhi]\:2606-\[CurlyPhi]p] (Z[z,zp]^2/(2 \[Rho]p^2)+ Log[Abs[Z[z,zp]]]) 
BdS\[CurlyPhi]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
(*On the azimuthal line*)
BdS\[Rho]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BdS\[Rho]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
BdS\[CurlyPhi]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BdS\[CurlyPhi]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BdSz1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := Cos[\[CurlyPhi]\:2606-\[CurlyPhi]] (ArcTanh[Sin[\[Phi][\[CurlyPhi],\[CurlyPhi]p]]]-2 Sin[\[Phi][\[CurlyPhi],\[CurlyPhi]p]]) Sign[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]-(2-Sqrt[2] Sqrt[1-Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]]) Sin[\[CurlyPhi]\:2606-\[CurlyPhi]] 
(*On the radial line*)
BdS\[Rho]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0
BdS\[Rho]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0 
BdS\[CurlyPhi]1[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0
BdS\[CurlyPhi]2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0
BdSz2[\[CurlyPhi]\:2606_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := -Sign[\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]]Log[Abs[\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]]] Sin[\[CurlyPhi]\:2606-\[CurlyPhi]p]


(* ::Section::Closed:: *)
(*Radial magnetisation*)


(*Rename and format result*)
CylB\[Rho][M_,P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[B\[Rho]Ana[M,P,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
B\[Rho]Ana[M_,P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, B\[Rho]AnaAxis[M,\[Rho]p,\[CurlyPhi]p,zp,z] , (M u0)/(4\[Pi]) Sum[(-1)^(m+n+q) B\[Rho]Summand[P,\[Rho]p[[m]],\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z] ,{m,1,2},{n,1,2},{q,1,2}]]
B\[Rho]AnaAxis[M_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (M u0)/(4\[Pi]) Sum[(-1)^(m+n+q) B\[Rho]SummandAxis[\[Rho]p[[m]],\[CurlyPhi]p[[q]],zp[[n]],z],{m,1,2},{n,1,2},{q,1,2}]

(*Geometry special cases*)
B\[Rho]Ana[M_,P_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, B\[Rho]AnaAxis[M,\[Rho]p,{0,2\[Pi]},zp,z],(M u0)/(4\[Pi]) Sum[(-1)^(m+n ) B\[Rho]SummandAS[P,\[Rho]p[[m]],\[Rho],\[CurlyPhi],zp[[n]],z],{m,1,2},{n,1,2}] ]
B\[Rho]AnaAxis[M_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (M u0)/(4\[Pi]) Sum[(-1)^(m+n) B\[Rho]SummandAxis[\[Rho]p[[m]],{0,2\[Pi]},zp[[n]],z],{m,1,2},{n,1,2}] 

(*Summands*)
B\[Rho]Summand[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {B\[Rho]S\[Rho]1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + B\[Rho]S\[Rho]2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], B\[Rho]S\[CurlyPhi]1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + B\[Rho]S\[CurlyPhi]2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], B\[Rho]Sz[P,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}

B\[Rho]S\[Rho]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := -(2/\[Rho]) ( \[Rho]p Z[z,zp] )/R[\[Rho],\[Rho]p,z,zp]  (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-(\[Rho] L[\[Rho],z,zp])/(Z[z,zp]^2 \[Rho]p) ( (S[\[Rho],\[Rho]p,z,zp]^2)/(a[\[Rho],z,zp]^2)  EllipticPiT[\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]+ (\!\(\*OverscriptBox[\(S\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2)/(\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2)  EllipticPiT[a[\[Rho],z,zp]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]))
B\[Rho]S\[Rho]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Which[(\[CurlyPhi]p==\[CurlyPhi]+\[Pi])|| (\[CurlyPhi]p==\[CurlyPhi]-\[Pi]),0,True,-Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]ArcTan[\[CapitalUpsilon][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]]-Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]ArcTanh[(G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1) /Z[z,zp]]]
B\[Rho]S\[CurlyPhi]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Z[z,zp] /\[Rho] Log[\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]+G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1]
B\[Rho]S\[CurlyPhi]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Which[(\[CurlyPhi]p==\[CurlyPhi]+\[Pi])|| (\[CurlyPhi]p==\[CurlyPhi]-\[Pi]),ArcTanh[R[\[Rho],\[Rho]p,z,zp]/Z[z,zp]],True,+ArcTan[\[CapitalUpsilon][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]] Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]-Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]] ArcTanh[(G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1) /Z[z,zp]]]
B\[Rho]Sz[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := -((2  \[Rho]p)/R[\[Rho],\[Rho]p,z,zp] EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]+\[CurlyPhi]p \[Alpha]1[\[Rho],\[Rho]p,z,zp]+ \[CurlyPhi]p \[Beta]1[\[Rho],\[Rho]p,z,zp,P] - \[Delta]1\[Zeta]\[Eta]\[Iota][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp,P])

B\[Rho]SummandAxis[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {B\[Rho]SAxisx[\[Rho]p,\[CurlyPhi]p,zp,z], B\[Rho]SAxisy[\[Rho]p,\[CurlyPhi]p,zp,z], B\[Rho]SAxisz[\[Rho]p,\[CurlyPhi]p,zp,z]}

B\[Rho]SAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] :=  Sin[\[CurlyPhi]p](Z[z,zp] /Sqrt[Z[z,zp] ^2+\[Rho]p^2]+ArcTanh[Sqrt[Z[z,zp]^2+\[Rho]p^2]/Z[z,zp] ])
B\[Rho]SAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -Cos[\[CurlyPhi]p](Z[z,zp]/Sqrt[Z[z,zp]^2+\[Rho]p^2]+ArcTanh[Sqrt[Z[z,zp]^2+\[Rho]p^2]/Z[z,zp] ])
B\[Rho]SAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := \[CurlyPhi]p(\[Rho]p /Sqrt[Z[z,zp]^2+\[Rho]p^2]- ArcTanh[\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2]])

B\[Rho]SummandAS[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := {B\[Rho]SAS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], B\[Rho]SAS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], B\[Rho]SASz[P,\[Rho]p,\[Rho],\[CurlyPhi],zp,z]}

B\[Rho]SAS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 4/\[Rho] ( \[Rho]p Z[z,zp])/R[\[Rho],\[Rho]p,z,zp]  (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]-(\[Rho] L[\[Rho],z,zp])/(Z[z,zp]^2 \[Rho]p) ( (S[\[Rho],\[Rho]p,z,zp]^2)/(a[\[Rho],z,zp]^2)  EllipticPi[\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2,k[\[Rho],\[Rho]p,z,zp]^2]+ (\!\(\*OverscriptBox[\(S\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2)/\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2 EllipticPi[a[\[Rho],z,zp]^2,k[\[Rho],\[Rho]p,z,zp]^2]))
B\[Rho]SAS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
B\[Rho]SASz[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := ((4\[Rho]p)/R[\[Rho],\[Rho]p,z,zp] EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]-2\[Pi] \[Alpha]1[\[Rho],\[Rho]p,z,zp]-2\[Pi] \[Beta]1[\[Rho],\[Rho]p,z,zp,P])
(*Along the axis*)
B\[Rho]SAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
B\[Rho]SAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
B\[Rho]SAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := -\[CurlyPhi]p Log[\[Rho]p]
(*Along the axis & axisymmetric*)
B\[Rho]SAxisx[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
B\[Rho]SAxisy[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
B\[Rho]SAxisz[\[Rho]p_,{0,2\[Pi]},zp_,z_] := B\[Rho]SAxisz[\[Rho]p,2\[Pi],zp,z]
B\[Rho]SAxisx[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
B\[Rho]SAxisy[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
B\[Rho]SAxisz[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := B\[Rho]SAxisz[\[Rho]p,2\[Pi],zp,zp] 
(*Solid*)
B\[Rho]S\[Rho]1[0,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := -(L[\[Rho],z,zp] /\[Rho])ArcTan[(\[Rho] Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]])/Z[z,zp]]
B\[Rho]S\[Rho]1[0,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
(*Solid & axisymmetric*)
B\[Rho]SAS\[Rho][0,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
(*On the shell plane*)
EllipticPiT[1,\[Phi]_,k_] := EllipticFT[\[Phi],k]-1/(1-k)(EllipticET[\[Phi],k]-Sqrt[1-k Sin[\[Phi]]^2] Tan[\[Phi]]);
(*On the shell plane & axisymmetric*)
B\[Rho]SAS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the section plane*)
B\[Rho]S\[Rho]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
B\[Rho]S\[CurlyPhi]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := -ArcTanh[\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]/Z[z,zp]]
(*On the disc plane*)
B\[Rho]S\[Rho]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
B\[Rho]S\[Rho]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
B\[Rho]S\[CurlyPhi]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
B\[Rho]S\[CurlyPhi]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the axial line*)
B\[Rho]S\[Rho]2[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
B\[Rho]S\[CurlyPhi]2[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := -Log[Abs[Z[z,zp]]] Sign[Z[z,zp]]
(*On the azimuthal line*)
B\[Rho]S\[Rho]1[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
B\[Rho]S\[Rho]2[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
B\[Rho]S\[CurlyPhi]1[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
B\[Rho]S\[CurlyPhi]2[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
EllipticFT[\[Phi]_,1]:=Sin[\[Phi]]CarlsonRC[1,Cos[\[Phi]]^2]
(*On the radial line*)
B\[Rho]S\[Rho]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0
B\[Rho]S\[Rho]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0 
B\[Rho]S\[CurlyPhi]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0
B\[Rho]S\[CurlyPhi]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0


(* ::Section::Closed:: *)
(*Azimuthal magnetisation*)


(*Rename and format result*)
CylB\[CurlyPhi][M_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[B\[CurlyPhi]Ana[M,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
B\[CurlyPhi]Ana[M_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, B\[CurlyPhi]AnaAxis[M,\[Rho]p,\[CurlyPhi]p,zp,z] , (M u0)/(4\[Pi]) (Sum[(-1)^(m+n+q) B\[CurlyPhi]Summand[\[Rho]p[[m]],\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z],{m,1,2},{n,1,2},{q,1,2}]+ M\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z])]
B\[CurlyPhi]AnaAxis[M_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (M u0)/(4\[Pi]) Sum[(-1)^(m+n+q) B\[CurlyPhi]SummandAxis[\[Rho]p[[m]],\[CurlyPhi]p[[q]],zp[[n]],z],{m,1,2},{n,1,2},{q,1,2}]

(*Magnetisation vector*)
M\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[InsideVolume[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],4\[Pi]{0,1,0},{0,0,0}]

(*Geometry special cases*)
B\[CurlyPhi]Ana[M_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0,B\[CurlyPhi]AnaAxis[M,\[Rho]p,{0,2\[Pi]},zp,z],(M u0)/(4\[Pi]) (M\[CurlyPhi][\[Rho]p,\[Rho],{0,2\[Pi]},\[Pi],zp,z])]
B\[CurlyPhi]AnaAxis[M_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := {0,0,0}

(*Summands*)
B\[CurlyPhi]Summand[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {B\[CurlyPhi]S\[Rho][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] , B\[CurlyPhi]S\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], B\[CurlyPhi]Sz[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}

B\[CurlyPhi]S\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Which[(\[CurlyPhi]p==\[CurlyPhi]+\[Pi]) || (\[CurlyPhi]p==\[CurlyPhi]-\[Pi]),ArcTanh[R[\[Rho],\[Rho]p,z,zp]/Z[z,zp]],True,ArcTan[\[CapitalUpsilon][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]]Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]-Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]] ArcTanh[G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1/Z[z,zp]]]
B\[CurlyPhi]S\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Which[(\[CurlyPhi]p==\[CurlyPhi]+\[Pi])|| (\[CurlyPhi]p==\[CurlyPhi]-\[Pi]),0,True,ArcTan[\[CapitalUpsilon][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]]Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]] + Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]ArcTanh[G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1/Z[z,zp]]]
B\[CurlyPhi]Sz[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Log[\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]+G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1]


B\[CurlyPhi]SummandAxis[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {B\[CurlyPhi]SAxisx[\[Rho]p,\[CurlyPhi]p,zp,z], B\[CurlyPhi]SAxisy[\[Rho]p,\[CurlyPhi]p,zp,z], B\[CurlyPhi]SAxisz[\[Rho]p,\[CurlyPhi]p,zp,z]}

B\[CurlyPhi]SAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -ArcTanh[Sqrt[Z[z,zp]^2+\[Rho]p^2]/Z[z,zp]]Cos[\[CurlyPhi]p]
B\[CurlyPhi]SAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -ArcTanh[Sqrt[Z[z,zp]^2+\[Rho]p^2]/Z[z,zp]]Sin[\[CurlyPhi]p]
B\[CurlyPhi]SAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := ArcTanh[\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2]]
(*Along the axis*)
B\[CurlyPhi]SAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
B\[CurlyPhi]SAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
B\[CurlyPhi]SAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := Log[\[Rho]p]
(*On the section plane*)
B\[CurlyPhi]S\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := -ArcTanh[\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]/Z[z,zp]]
B\[CurlyPhi]S\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
(*On the disc plane*)
B\[CurlyPhi]S\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
B\[CurlyPhi]S\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the axial line*)
B\[CurlyPhi]S\[Rho][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := -Log[Abs[Z[z,zp]]]Sign[Z[z,zp]]
B\[CurlyPhi]S\[CurlyPhi][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
(*On the azimuthal line*)
B\[CurlyPhi]S\[Rho][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
B\[CurlyPhi]S\[CurlyPhi][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the radial line*)
B\[CurlyPhi]S\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0 
B\[CurlyPhi]S\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0
B\[CurlyPhi]Sz[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := -Sign[\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]]Log[Abs[\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]]] 



(* ::Section::Closed:: *)
(*Axial magnetisation*)


(*Rename and format result*)
CylB\[ScriptZ][M_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[BzAna[M,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
BzAna[M_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BzAnaAxis[M,\[Rho]p,\[CurlyPhi]p,zp,z] , (M u0)/(4\[Pi]) Sum[(-1)^(m+n+q) BzSummand[\[Rho]p[[m]],\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z] ,{m,1,2},{n,1,2},{q,1,2}]]
BzAnaAxis[M_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (M u0)/(4\[Pi]) Sum[(-1)^(m+n+q) BzSummandAxis[\[Rho]p[[m]],\[CurlyPhi]p[[q]],zp[[n]],z],{m,1,2},{n,1,2},{q,1,2}]

(*Geometry special cases*)
BzAna[M_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BzAnaAxis[M,\[Rho]p,{0,2\[Pi]},zp,z],(M u0)/(4\[Pi]) Sum[(-1)^(m+n ) BzSummandAS[\[Rho]p[[m]],\[Rho],\[CurlyPhi],zp[[n]],z],{m,1,2},{n,1,2}] ]
BzAnaAxis[M_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (M u0)/(4\[Pi]) Sum[(-1)^(m+n) BzSummandAxis[\[Rho]p[[m]],{0,2\[Pi]},zp[[n]],z],{m,1,2},{n,1,2}] 
BzAna[M_,{0,\[Rho]p_},\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := (M u0)/(4\[Pi]) Sum[(-1)^(q+n) (BzSummand1[\[Rho]p,\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z] + Sum[(-1)^m BzSummand2[(m-1)\[Rho]p,\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z],{m,1,2}]),{n,1,2},{q,1,2}]
BzAna[M_,{0,\[Rho]p_},\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := (M u0)/(4\[Pi]) Sum[(-1)^n BzSummandAS[\[Rho]p,\[Rho],\[CurlyPhi],zp[[n]],z] ,{n,1,2}]

(*Summands*)
BzSummand[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BzS\[Rho]1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + BzS\[Rho]2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BzS\[CurlyPhi]1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + BzS\[CurlyPhi]2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BzSz1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z] + BzSz2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}
BzSummand1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BzS\[Rho]1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BzS\[CurlyPhi]1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BzSz1[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}
BzSummand2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BzS\[Rho]2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BzS\[CurlyPhi]2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BzSz2[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}

BzS\[Rho]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := (2 \[Rho]p)/R[\[Rho],\[Rho]p,zp,z] (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-2EllipticDT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2])
BzS\[Rho]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Log[\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]+G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1]Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]
BzS\[CurlyPhi]1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := 1/\[Rho]  G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1 
BzS\[CurlyPhi]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Log[\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]+G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1]Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]
BzSz1[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Z[z,zp]/R[\[Rho],\[Rho]p,zp,z] (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]+ (\[Rho]p^2-\[Rho]^2)/\[CurlyRho][\[Rho],\[Rho]p]^2 EllipticPiT[\[Kappa][\[Rho],\[Rho]p]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]) 
BzSz2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Which[(\[CurlyPhi]p==\[CurlyPhi]+\[Pi])|| (\[CurlyPhi]p==\[CurlyPhi]-\[Pi]),0,True,-ArcTan[\[CapitalUpsilon][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]]]

BzSummandAxis[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {BzSAxisx[\[Rho]p,\[CurlyPhi]p,zp,z], BzSAxisy[\[Rho]p,\[CurlyPhi]p,zp,z], BzSAxisz[\[Rho]p,\[CurlyPhi]p,zp,z]}

BzSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] :=  Sin[\[CurlyPhi]p](\[Rho]p /Sqrt[Z[z,zp]^2+\[Rho]p^2]-ArcTanh[\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2]] )
BzSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -Cos[\[CurlyPhi]p](\[Rho]p /Sqrt[Z[z,zp]^2+\[Rho]p^2]-ArcTanh[\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2]] )
BzSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := - ((Z[z,zp] \[CurlyPhi]p)/Sqrt[Z[z,zp]^2+\[Rho]p^2])


BzSummandAS[\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := {BzSAS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BzSAS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BzSASz[\[Rho]p,\[Rho],\[CurlyPhi],zp,z]}

BzSAS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := (-4  \[Rho]p)/R[\[Rho],\[Rho]p,zp,z]  (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]-2EllipticD[k[\[Rho],\[Rho]p,z,zp]^2])
BzSAS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
BzSASz[\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := (-2 Z[z,zp])/R[\[Rho],\[Rho]p,zp,z] (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]+ (\[Rho]p^2-\[Rho]^2)/\[CurlyRho][\[Rho],\[Rho]p]^2 EllipticPi[\[Kappa][\[Rho],\[Rho]p]^2,k[\[Rho],\[Rho]p,z,zp]^2])
(*Along the axis*)
BzSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := Sin[\[CurlyPhi]p](1-Log[\[Rho]p] )
BzSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := Cos[\[CurlyPhi]p](Log[\[Rho]p]-1)
BzSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
(*Along the axis & axisymmetric*)
BzSAxisx[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BzSAxisy[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BzSAxisz[\[Rho]p_,{0,2\[Pi]},zp_,z_] := BzSAxisz[\[Rho]p,2\[Pi],zp,z]
BzSAxisx[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
BzSAxisy[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
BzSAxisz[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
(*On the shell plane*)
BzSz1[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Z[z,zp] /Sqrt[Z[z,zp]^2+4 \[Rho]p^2]  EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],(4 \[Rho]p^2)/(Z[z,zp]^2+4 \[Rho]p^2)]
(*On the shell plane & axisymmetric *)
BzSASz[\[Rho]p_,\[Rho]p_,\[CurlyPhi]_,zp_,z_] := (-2 Z[z,zp])/R[\[Rho]p,\[Rho]p,zp,z] EllipticK[k[\[Rho]p,\[Rho]p,z,zp]^2]
(*On the section plane*)
BzSz2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
(*On the axial line*)
BzSz2[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,z_] := 0
(*On the azimuthal line*)
EllipticFT[\[Phi]_,1]:=Sin[\[Phi]]CarlsonRC[1,Cos[\[Phi]]^2]
BzSz1[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the radial line*)
BzS\[Rho]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0 
BzS\[CurlyPhi]2[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := -Sign[\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]]Log[Abs[\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]]] 


(* ::Section::Closed:: *)
(*Filament with azimuthal current density*)


(*Rename and format result*)
CylB\[ScriptCapitalI]\[ScriptF][I_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[BiAna[I,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
BiAna[I_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BiAnaAxis[I,\[Rho]p,\[CurlyPhi]p,zp,z] , (I u0)/(4\[Pi]) Sum[(-1)^q BiSummand[\[Rho]p,\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp,z],{q,1,2}]]
BiAnaAxis[I_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (I u0)/(4\[Pi]) Sum[(-1)^q BiSummandAxis[\[Rho]p,\[CurlyPhi]p[[q]],zp,z],{q,1,2}]

(*Geometry special cases*)
BiAna[I_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0,BiAnaAxis[I,\[Rho]p,{0,2\[Pi]},zp,z],(I u0)/(4\[Pi]) BiSummandAS[\[Rho]p,\[Rho],\[CurlyPhi],zp,z]]
BiAnaAxis[I_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (I u0)/(4\[Pi]) BiSummandAxis[\[Rho]p,{0,2\[Pi]},zp,z]

(*Summands*)
BiSummand[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BiS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BiS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BiSz[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}

BiS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Z[z,zp]/\[Rho]  (1/R[\[Rho],\[Rho]p,z,zp] (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-(1+1/2 \!\(\*OverscriptBox[\(k\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2)EllipticET[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho], \[Rho]p,z,zp]^2])-1/2 (k[\[Rho],\[Rho]p,z,zp]^2 T[\[Rho],\[Rho]p,z,zp]^2)/\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2  Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]/ (G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1))
BiS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := -(Z[z,zp]/\[Rho] )G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp] 
BiSz[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := -(1/R[\[Rho], \[Rho]p,z,zp])(EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho], \[Rho]p,z,zp]^2]-\!\(\*OverscriptBox[\(T\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2/\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2 EllipticET[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho], \[Rho]p,z,zp]^2])+1/2 (k[\[Rho],\[Rho]p,z,zp]^2 \!\(\*OverscriptBox[\(T\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2)/\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2  Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]/ (G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1)

BiSummandAxis[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {BiSAxisx[\[Rho]p,\[CurlyPhi]p,zp,z], BiSAxisy[\[Rho]p,\[CurlyPhi]p,zp,z], BiSAxisz[\[Rho]p,\[CurlyPhi]p,zp,z]}

BiSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (Z[z,zp] \[Rho]p Sin[\[CurlyPhi]p])/(Z[z,zp]^2+\[Rho]p^2)^(3/2)
BiSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -((Z[z,zp] \[Rho]p Cos[\[CurlyPhi]p])/(Z[z,zp]^2+\[Rho]p^2)^(3/2))
BiSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (\[Rho]p^2 \[CurlyPhi]p)/(Z[z,zp]^2+\[Rho]p^2)^(3/2)

BiSummandAS[\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := {BiSAS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BiSAS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BiSASz[\[Rho]p,\[Rho],\[CurlyPhi],zp,z]}

BiSAS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := (2 Z[z,zp])/(\[Rho] R[\[Rho]p,\[Rho],z,zp]) (T[\[Rho],\[Rho]p,z,zp]^2/\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2 EllipticE[k[\[Rho]p,\[Rho],z,zp]^2]-EllipticK[k[\[Rho]p,\[Rho],z,zp]^2]) 
BiSAS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
BiSASz[\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 2/R[\[Rho],\[Rho]p,z,zp] (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]-\!\(\*OverscriptBox[\(T\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2/\!\(\*OverscriptBox[\(R\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2 EllipticE[k[\[Rho],\[Rho]p,z,zp]^2] )

(*Along the axis & axisymmetric*)
BiSAxisx[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BiSAxisy[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BiSAxisz[\[Rho]p_,{0,2\[Pi]},zp_,z_] := BiSAxisz[\[Rho]p,2\[Pi],zp,z]
(*On the azimuthal line*)
BiS\[Rho][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BiS\[CurlyPhi][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BiSz[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := -Sign[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]] ArcTanh[Sin[\[Phi][\[CurlyPhi],\[CurlyPhi]p]]] /(2 \[Rho]p)


(* ::Section::Closed:: *)
(*Disc with azimuthal current density*)


(*Rename and format result*)
CylB\[ScriptCapitalK]\[ScriptD][K_,P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[BkAna[K,P,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
BkAna[K_,P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BkAnaAxis[K,\[Rho]p,\[CurlyPhi]p,zp,z] , (K u0)/(4\[Pi]) Sum[(-1)^(m+q) BkSummand[P,\[Rho]p[[m]],\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp,z],{m,1,2},{q,1,2}]]
BkAnaAxis[K_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (K u0)/(4\[Pi]) Sum[(-1)^(m+q) BkSummandAxis[\[Rho]p[[m]],\[CurlyPhi]p[[q]],zp,z],{m,1,2},{q,1,2}]

(*Geometry special cases*)
BkAna[K_,P_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] :=If[\[Rho]==0,BkAnaAxis[K,\[Rho]p,{0,2\[Pi]},zp,z], (K u0)/(4\[Pi]) Sum[(-1)^m BkSummandAS[P,\[Rho]p[[m]],\[Rho],\[CurlyPhi],zp,z] ,{m,1,2}]]
BkAnaAxis[K_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (K u0)/(4\[Pi]) Sum[(-1)^m BkSummandAxis[\[Rho]p[[m]],{0,2\[Pi]},zp,z],{m,1,2}]

(*Summands*)
BkSummand[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BkS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BkS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BkSz[P,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}

BkS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := 2/\[Rho] ( \[Rho]p Z[z,zp] )/R[\[Rho],\[Rho]p,z,zp]  (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-(\[Rho] L[\[Rho],z,zp])/(Z[z,zp]^2 \[Rho]p) ( (S[\[Rho],\[Rho]p,z,zp]^2)/(a[\[Rho],z,zp]^2)  EllipticPiT[\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]+ (\!\(\*OverscriptBox[\(S\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2)/(\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2)  EllipticPiT[a[\[Rho],z,zp]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]))
BkS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := -(Z[z,zp] /\[Rho])Log[\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]+G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1]
BkSz[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := (2  \[Rho]p)/R[\[Rho],\[Rho]p,z,zp] EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]+\[CurlyPhi]p \[Alpha]1[\[Rho],\[Rho]p,z,zp]+ \[CurlyPhi]p \[Beta]1[\[Rho],\[Rho]p,z,zp,P] - \[Delta]1\[Zeta]\[Eta]\[Iota][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp,P]

BkSummandAxis[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {BkSAxisx[\[Rho]p,\[CurlyPhi]p,zp,z], BkSAxisy[\[Rho]p,\[CurlyPhi]p,zp,z], BkSAxisz[\[Rho]p,\[CurlyPhi]p,zp,z]}

BkSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -Sin[\[CurlyPhi]p] Z[z,zp] /Sqrt[Z[z,zp] ^2+\[Rho]p^2]
BkSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := Cos[\[CurlyPhi]p] Z[z,zp]/Sqrt[Z[z,zp]^2+\[Rho]p^2]
BkSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -\[CurlyPhi]p(\[Rho]p /Sqrt[Z[z,zp]^2+\[Rho]p^2]- ArcTanh[\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2]])

BkSummandAS[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := {BkSAS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BkSAS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BkSASz[P,\[Rho]p,\[Rho],\[CurlyPhi],zp,z]}

BkSAS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := -(4/\[Rho]) ( \[Rho]p Z[z,zp])/R[\[Rho],\[Rho]p,z,zp]  (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]-(\[Rho] L[\[Rho],z,zp])/(Z[z,zp]^2 \[Rho]p) ( (S[\[Rho],\[Rho]p,z,zp]^2)/(a[\[Rho],z,zp]^2)  EllipticPi[\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2,k[\[Rho],\[Rho]p,z,zp]^2]+ (\!\(\*OverscriptBox[\(S\), \(_\)]\)[\[Rho],\[Rho]p,z,zp]^2)/\!\(\*OverscriptBox[\(a\), \(_\)]\)[\[Rho],z,zp]^2 EllipticPi[a[\[Rho],z,zp]^2,k[\[Rho],\[Rho]p,z,zp]^2]))
BkSAS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
BkSASz[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := -((4\[Rho]p)/R[\[Rho],\[Rho]p,z,zp])EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]+2\[Pi] \[Alpha]1[\[Rho],\[Rho]p,z,zp]+2\[Pi] \[Beta]1[\[Rho],\[Rho]p,z,zp,P]

(*Along the axis*)
BkSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
BkSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
BkSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := \[CurlyPhi]p Log[\[Rho]p]
(*Along the axis & axisymmetric*)
BkSAxisx[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BkSAxisy[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BkSAxisz[\[Rho]p_,{0,2\[Pi]},zp_,z_] := BkSAxisz[\[Rho]p,2\[Pi],zp,z]
(*Solid*)
BkS\[Rho][0,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := L[\[Rho],z,zp] /\[Rho] ArcTan[(\[Rho] Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]])/Z[z,zp]]
BkS\[Rho][0,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
(*Solid & axisymmetric*)
BkSAS\[Rho][0,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
(*On the disc plane*)
BkS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BkS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0 
(*On the disc plane & axisymmetric*)
BkSAS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the azimuthal line*)
EllipticFT[\[Phi]_,1]:=Sin[\[Phi]]CarlsonRC[1,Cos[\[Phi]]^2]
BkS\[Rho][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
BkS\[CurlyPhi][\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,zp_] := 0
(*On the radial line*)
BkS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0
BkS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := 0 


(* ::Section::Closed:: *)
(*Shell with azimuthal current density*)


(*Rename and format result*)
CylB\[ScriptCapitalK]\[ScriptS][K_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[BsAna[K,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
BsAna[K_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BsAnaAxis[K,\[Rho]p,\[CurlyPhi]p,zp,z] , (K u0)/(4\[Pi]) Sum[(-1)^(n+q) BsSummand[\[Rho]p,\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z] ,{n,1,2},{q,1,2}]]
BsAnaAxis[K_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (K u0)/(4\[Pi]) Sum[(-1)^(n+q) BsSummandAxis[\[Rho]p,\[CurlyPhi]p[[q]],zp[[n]],z],{n,1,2},{q,1,2}]

(*Geometry special cases*)
BsAna[K_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_]:= If[\[Rho]==0, BsAnaAxis[K,\[Rho]p,{0,2\[Pi]},zp,z],(K u0)/(4\[Pi]) Sum[(-1)^n  BsSummandAS[\[Rho]p,\[Rho],\[CurlyPhi],zp[[n]],z],{n,1,2}] ]
BsAnaAxis[K_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (K u0)/(4\[Pi]) Sum[(-1)^n BsSummandAxis[\[Rho]p,{0,2\[Pi]},zp[[n]],z],{n,1,2}] 

(*Summands*)
BsSummand[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BsS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BsS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BsSz[\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]}

BsS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := (2 \[Rho]p)/R[\[Rho],\[Rho]p,zp,z] (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]-2EllipticDT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2])
BsS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := 1/\[Rho]  G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1 
BsSz[\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Z[z,zp]/R[\[Rho],\[Rho]p,zp,z] (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]+ (\[Rho]p^2-\[Rho]^2)/\[CurlyRho][\[Rho],\[Rho]p]^2 EllipticPiT[\[Kappa][\[Rho],\[Rho]p]^2,\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2])

BsSummandAxis[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {BsSAxisx[\[Rho]p,\[CurlyPhi]p,zp,z], BsSAxisy[\[Rho]p,\[CurlyPhi]p,zp,z], BsSAxisz[\[Rho]p,\[CurlyPhi]p,zp,z]}

BsSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] :=  Sin[\[CurlyPhi]p] \[Rho]p /Sqrt[Z[z,zp]^2+\[Rho]p^2]
BsSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -Cos[\[CurlyPhi]p] \[Rho]p /Sqrt[Z[z,zp]^2+\[Rho]p^2]
BsSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := - ((Z[z,zp] \[CurlyPhi]p)/Sqrt[Z[z,zp]^2+\[Rho]p^2])

BsSummandAS[\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := {BsSAS\[Rho][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BsSAS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BsSASz[\[Rho]p,\[Rho],\[CurlyPhi],zp,z]}

BsSAS\[Rho][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := (4  \[Rho]p)/R[\[Rho],\[Rho]p,zp,z]  (2EllipticD[k[\[Rho],\[Rho]p,z,zp]^2]-EllipticK[k[\[Rho],\[Rho]p,z,zp]^2])
BsSAS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
BsSASz[\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := (-2 Z[z,zp])/R[\[Rho],\[Rho]p,zp,z] (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]+ (\[Rho]p^2-\[Rho]^2)/\[CurlyRho][\[Rho],\[Rho]p]^2 EllipticPi[\[Kappa][\[Rho],\[Rho]p]^2,k[\[Rho],\[Rho]p,z,zp]^2])

(*Along the axis*)
BsSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := Sin[\[CurlyPhi]p]
BsSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := -Cos[\[CurlyPhi]p]
BsSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
(*Along the axis & axisymmetric*)
BsSAxisx[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BsSAxisy[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BsSAxisz[\[Rho]p_,{0,2\[Pi]},zp_,z_] := BsSAxisz[\[Rho]p,2\[Pi],zp,z]
BsSAxisx[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
BsSAxisy[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
BsSAxisz[\[Rho]p_,{0,2\[Pi]},zp_,zp_] := 0
(*On the shell plane*)
BsSz[\[Rho]p_,\[Rho]p_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := Z[z,zp] /Sqrt[Z[z,zp]^2+4 \[Rho]p^2]  EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],(4 \[Rho]p^2)/(Z[z,zp]^2+4 \[Rho]p^2)]
(*On the shell plane & axisymmetric*)
BsSASz[\[Rho]p_,\[Rho]p_,\[CurlyPhi]_,zp_,z_] := (-2 Z[z,zp])/R[\[Rho]p,\[Rho]p,zp,z] EllipticK[k[\[Rho]p,\[Rho]p,z,zp]^2]
(*On the azimuthal line*)
EllipticFT[\[Phi]_,1]:=Sin[\[Phi]]CarlsonRC[1,Cos[\[Phi]]^2]


(* ::Section::Closed:: *)
(*Volume with azimuthal current density*)


(*Rename and format result*)
CylB\[ScriptCapitalJ]\[ScriptV][J_,P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_]:=Chop[N[BcAna[J,P,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z],$MachinePrecision]]

(*Function handles*)
BcAna[J_,P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BcAnaAxis[J,\[Rho]p,\[CurlyPhi]p,zp,z], (J u0)/(4\[Pi]) Sum[(-1)^(m+n+q) BcSummand[P,\[Rho]p[[m]],\[Rho],\[CurlyPhi]p[[q]],\[CurlyPhi],zp[[n]],z] ,{m,1,2},{n,1,2},{q,1,2}]]
BcAnaAxis[J_,\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := (J u0)/(4\[Pi]) Sum[(-1)^(m+n+q) BcSummandAxis[\[Rho]p[[m]],\[CurlyPhi]p[[q]],zp[[n]],z],{m,1,2},{n,1,2},{q,1,2}]

(*Geometry special cases*)
BcAna[J_,P_,\[Rho]p_,\[Rho]_,{0,2\[Pi]},\[CurlyPhi]_,zp_,z_] := If[\[Rho]==0, BcAnaAxis[J,\[Rho]p,{0,2\[Pi]},zp,z],(J u0)/(4\[Pi]) Sum[(-1)^(m+n ) BcSummandAS[P,\[Rho]p[[m]],\[Rho],\[CurlyPhi],zp[[n]],z],{m,1,2},{n,1,2}] ]
BcAnaAxis[J_,\[Rho]p_,{0,2\[Pi]},zp_,z_] := (J u0)/(4\[Pi]) Sum[(-1)^(m+n) BcSummandAxis[\[Rho]p[[m]],{0,2\[Pi]},zp[[n]],z],{m,1,2},{n,1,2}] 

(*Summands*)
BcSummand[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := {BcS\[Rho][P,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BcS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z], BcSz[P,\[Rho]p,\[Rho],\[CurlyPhi]p,\[CurlyPhi],zp,z]} 

BcS\[Rho][P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := 1/2 (R[\[Rho],\[Rho]p,z,zp] 4/3 (EllipticFT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2]+(-2+k[\[Rho],\[Rho]p,z,zp]^2)EllipticDT[\[Phi][\[CurlyPhi],\[CurlyPhi]p],k[\[Rho],\[Rho]p,z,zp]^2])-4/3 Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1-\[Rho] (\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]+1/2 Sin[2 \[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]])( \[Alpha]1[\[Rho],\[Rho]p,z,zp] + \[Beta]1[\[Rho],\[Rho]p,z,zp,P])-\[Rho] \[Delta]3[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp,P])
BcS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := 1/(2\[Rho] ) ( (\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]) G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1+ (Z[z,zp]^2+\[Rho]^2 Sin[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]^2) Log[\[Rho]p-\[Rho] Cos[\[CapitalPhi][\[CurlyPhi],\[CurlyPhi]p]]+G[\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp]^-1])
BcSz[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]_,zp_,z_] := \[CurlyPhi]p \[Alpha]3[\[Rho],\[Rho]p,z,zp]+1/2 Sign[Z[z,zp]]\[Rho]p(2\[CurlyPhi]p \[Alpha]2[\[Rho],\[Rho]p,z,zp]+\[CurlyPhi]p \[Beta]2[\[Rho],\[Rho]p,z,zp,P]-\[CurlyPhi]p Sqrt[\[Pi]]/8 \[Beta]3[\[Rho],\[Rho]p,z,zp,P] + \[Gamma]3\[Zeta]\[Eta]\[Iota][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp,P] - \[Delta]2\[Zeta]\[Eta]\[Iota][\[Rho],\[Rho]p,\[CurlyPhi],\[CurlyPhi]p,z,zp,P])

BcSummandAxis[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := {BcSAxisx[\[Rho]p,\[CurlyPhi]p,zp,z], BcSAxisy[\[Rho]p,\[CurlyPhi]p,zp,z], BcSAxisz[\[Rho]p,\[CurlyPhi]p,zp,z]}

BcSAxisx[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := Sin[\[CurlyPhi]p]Sqrt[Z[z,zp] ^2+\[Rho]p^2]
BcSAxisy[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -Cos[\[CurlyPhi]p] Sqrt[Z[z,zp] ^2+\[Rho]p^2]
BcSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,z_] := -\[CurlyPhi]p Z[z,zp] ArcTanh[\[Rho]p/Sqrt[Z[z,zp]^2+\[Rho]p^2]]

BcSummandAS[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := {BcSAS\[Rho][P,\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BcSAS\[CurlyPhi][\[Rho]p,\[Rho],\[CurlyPhi],zp,z], BcSASz[P,\[Rho]p,\[Rho],\[CurlyPhi],zp,z]}

BcSAS\[Rho][P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := \[Pi] \[Rho] (\[Alpha]1[\[Rho],\[Rho]p,z,zp]+\[Beta]1[\[Rho],\[Rho]p,z,zp,P]+ \[CurlyTheta][\[Rho],\[Rho]p,z,zp,P])- R[\[Rho],\[Rho]p,z,zp] 4/3 (EllipticK[k[\[Rho],\[Rho]p,z,zp]^2]+(-2+k[\[Rho], \[Rho]p,z,zp]^2)EllipticD[k[\[Rho], \[Rho]p,z,zp]^2])
BcSAS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 0
BcSASz[P_,\[Rho]p_,\[Rho]_,\[CurlyPhi]_,zp_,z_] := 2\[Pi](\[Alpha]3[\[Rho],\[Rho]p,z,zp]+1/2 \[Rho]p Sign[Z[z,zp]](2 \[Alpha]2[\[Rho],\[Rho]p,z,zp]+\[Beta]2[\[Rho],\[Rho]p,z,zp,P]-Sqrt[\[Pi]]/8 \[Beta]3[\[Rho],\[Rho]p,z,zp,P]))
(*Along the axis*)
BcSAxisz[\[Rho]p_,\[CurlyPhi]p_,zp_,zp_] := 0
(*Along the axis & axisymmetric*)
BcSAxisx[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BcSAxisy[\[Rho]p_,{0,2\[Pi]},zp_,z_] := 0
BcSAxisz[\[Rho]p_,{0,2\[Pi]},zp_,z_] := BcSAxisz[\[Rho]p,2\[Pi],zp,z]
(*On the azimuthal line*)
EllipticFT[\[Phi]_,1]:=Sin[\[Phi]]CarlsonRC[1,Cos[\[Phi]]^2]
(*On the radial line*)
BcS\[CurlyPhi][\[Rho]p_,\[Rho]_,\[CurlyPhi]p_,\[CurlyPhi]p_,zp_,zp_] := (*-((Abs[Overscript[\[CurlyRho], _][\[Rho],\[Rho]p]] Overscript[\[CurlyRho], _][\[Rho],\[Rho]p])/(2 \[Rho]))*)-((Abs[\!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]] \[Rho]p (\[Rho]p-2 \[Rho]))/(2 \[Rho] \!\(\*OverscriptBox[\(\[CurlyRho]\), \(_\)]\)[\[Rho],\[Rho]p]))


(* ::Section:: *)
(**)


End[];


EndPackage[];
