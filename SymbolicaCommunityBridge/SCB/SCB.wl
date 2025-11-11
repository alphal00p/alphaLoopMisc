(* ::Package:: *)

(* ===== fast reload trap: if re-Get happens with SCB$Reloading=True, just reattach impls and exit ===== *)

If[ TrueQ[Global`SCB$Reloading] && NameQ["SCB`Private`iLoad"],
  SCB`InfoMessage["Reloading SymbolicaCommunityBridge function definitions..."];
  (* 
  ToExpression["SCB`Private`iLoad", StandardForm, Function[f, f["reload"]]];
  Global`SCB$Reloading = False;
  Throw[Null, "SCB/Reload"];
  *)
];

(* ===== public section ===== *)
BeginPackage["SCB`", {"GeneralUtilities`"}];

(* exported symbol registry *)
SCB`$PublicSymbols::usage = "List of exported SCB symbols.";
SCB`$GatedSymbols::usage  = "Subset of exported symbols that require Init[].";

Init::usage          = "SCB`Init[...] must be called before use.";
Reset::usage         = "SCB`Reset[] clears initialization.";
InitializedQ::usage  = "SCB`InitializedQ[] gives True if initialized.";
SimplifyColor::usage = "SCB`SimplifyColor[expr] simplifies color algebra.";
SimplifyGamma::usage = "SCB`SimplifyGamma[expr] simplifies gamma algebra.";
GammaLoop::usage     = "SCB`GammaLoop[expr, opts] executes commands and returns results.";
LoadSCB::usage       = "SCB`LoadSCB[] hot-reloads private implementations without restarting.";

Pkg::notinit = "`1` requires `2`.";
Init::req    = "Missing required options: `1`.";
Cmd::error        = "Error in SymbolicaCommunityBridge: `1`.";

Info::info = "`1`";
(* message text, user can disable with: Off[Pkg::welcome] *)
Welcome::info = "

Symbolica Community Bridge:
---------------------

Requirements:

* Python3, with: 
    
    'python3 -m pip install pydot'
    'python3 -m pip install pyzmq'

* Symbolica community: (last tested with revision: 5e40c2c2e3cfd5d60f2630524350a783582e8daa)
    
    'python3 -m pip install symbolica'
    
    or
    
    'git clone git@github.com:alphal00p/symbolica-community.git && cd symbolica-community'
    'cargo run --features \"python_stubgen\"'
    'maturin build --release --features \"python_stubgen\"'
    'python3 -m pip install <path_to_the_wheel_generated_above>

* GammaLoop API:

    'git clone git@github.com:alphal00p/gammaloop.git && cd gammaloop'
    'maturin develop -m gammaloop-api/Cargo.toml --features=ufo_support,python_api'

* Then, initialize SymbolicaCommunityBridge with:

SCB`.`Init[
 PythonInterpreter -> \".../.venv/bin/python3\", (* Or Automatic to rely on system defaults. *)                               
 GammaLoopPath -> \".../gammaloop_hedge_numerator/gammaloop-api/python\", (* Or Null to disable gammaLoop features *)
 SymbolicaCommunityPath -> \".../.venv/lib/python3.12/site-packages\", (* Or Null to disable symbolica community features *)
 SymbolicaLicense -> \"YourSymbolicLicense\" (* Optional *)
 FORMPath -> \".../bin/form\", (* Or Null to disable symbolica FORM features *)
 DebugLevel -> Automatic (* Or any string default to pass as tracing directive *)
];
";

(* declare once; add new API names here only *)
SCB`$PublicSymbols = {Init, Reset, InitializedQ, SimplifyColor, SimplifyGamma, GammaLoop, RunPython, ResetGammaLoop,SExpand,SetDebugLevel};
SCB`$GatedSymbols  = {SimplifyColor, SimplifyGamma, GammaLoop, RunPython,SExpand};

(* ===== private section ===== *)
Begin["`Private`"];

(* remember this file\[CloseCurlyQuote]s path for self-reload *)
If[!ValueQ[SCB`$ThisFile], SCB`$ThisFile = $InputFileName];

(* persistent state: initialize only if not yet set *)
If[!ValueQ[$initialized], $initialized = False];
If[!ValueQ[$config],      $config      = <||>];

(* realized symbol lists *)
publicSyms[] := SCB`$PublicSymbols;
gatedSyms[]  := SCB`$GatedSymbols;

(* protect/unprotect helpers for all exports *)
iUnprotectPublic[] := Quiet @ Unprotect @@ publicSyms[]; 
iProtectPublic[]   := Quiet @ Protect   @@ publicSyms[];

(* clear only definitions; keep usage/options/messages intact *)
iClearImpl[] := Scan[
  Function[s, (DownValues[s]=.; SubValues[s]=.; UpValues[s]=.; OwnValues[s]=.; FormatValues[s]=.)],
  publicSyms[]
];

SetAttributes[requireInit, HoldFirst];

iDefineImpl[] := Module[{},


  Options[Init]      = {PythonInterpreter -> Automatic, GammaLoopPath -> Automatic, GammaLoopStatePath -> Automatic, SymbolicaCommunityPath -> Automatic, FORMPath -> Automatic, SymbolicaLicense -> Automatic, DebugLevel -> Automatic, ParsingModeFullForm -> True};
  Options[GammaLoop] = {DebugLevel -> "Info"};

  (* guard for pre-init calls *)
  requireInit[s_Symbol] :=
    ( s[args___] /; !TrueQ[$initialized] :=
        ( Message[
            Pkg::notinit,
            HoldForm[s],
            HoldForm @ SCB`Init[GammaLoopPath -> "...", SymbolicaCommunityPath -> "..."]
          ];
          Null
        )
    );
  Scan[requireInit, gatedSyms[]];

  SCB`WrappedExternalEvaluate[py_,cmd_]:=Block[{res},
    If[$config["DebugLevel"]>0,
      InfoMessage["Running the following Python code:\n"<>cmd];
    ];
    res = ExternalEvaluate[py,cmd];
    If[Head[res]===Failure,
       Print[res];
       Throw[Null, "PythonError"];
    ];
    res
  ];
  
  (*
  SCB`UncurryExpr[expr_] := Module[{}, expr //. {x_[y___][z___] :> x[y, SCBridge`CurryDivider, z]} ] ;
  
  SCB`CurryExpr[expr_] := Module[{}, expr //. {x_[y___, SCBridge`CurryDivider, z___] :> x[y][z]} ];
  
  SetAttributes[SCB`UncurryExprHeld, HoldAllComplete];
  SCB`UncurryExprHeld[expr_] :=
  ReleaseHold @ FixedPoint[
    Replace[#, x_[y___][z___] :> x[y, SCBridge`CurryDivider, z], {0, Infinity}, Heads -> True] &,
    HoldComplete[expr]
  ];

  SetAttributes[SCB`CurryExprHeld, HoldAllComplete];
  SCB`CurryExprHeld[expr_] :=
  ReleaseHold @ FixedPoint[
    Replace[#, x_[y___, SCBridge`CurryDivider, z___] :> x[y][z], {0, Infinity}, Heads -> True] &,
    HoldComplete[expr]
  ];
  *)
  
  SCB`ExprToString[expr_, OptionsPattern[{BackTickReplace -> False, FullFormParsing -> Automatic}]]:=Block[{res, fullFormParsing},
     fullFormParsing = If[OptionValue[FullFormParsing] === Automatic,
        $config["ParsingModeFullForm"],
        OptionValue[FullFormParsing]
     ];
     res = If[StringQ[expr],
       expr,
       If[fullFormParsing,
         ToString[FullForm[expr, NumberMarks -> False]],
         ToString[InputForm[expr]]
       ]
     ];
     If[
     OptionValue[BackTickReplace],
     StringReplace[res,{"`"->"MATHEMATICABACKTICK"}],
     res
     ]
  ];

  SCB`StringToExpr[expr_, OptionsPattern[{BackTickReplace -> False}]]:=Block[{processedExpr,res},
  
     If[StringQ[expr],
       processedExpr = If[OptionValue[BackTickReplace],
          StringReplace[expr,{"MATHEMATICABACKTICK"->"`"}]
          ,
          expr
       ];
       res=If[And[StringLength[processedExpr]>=2, StringTake[processedExpr,1]==="'", StringTake[processedExpr,-1]==="'"],
         ToExpression[StringTake[processedExpr,{2,-2}]]
         ,
         ToExpression[processedExpr]
       ];
       res
       ,
       processedExpr
     ]
  ];

  SCB`InfoMessage[arg_] := Module[{},
	Message[SCB`Info::info, arg];
	Unprotect[$MessageList];
    $MessageList = DeleteCases[$MessageList, HoldForm[SCB`Info::info]];
    Protect[$MessageList];
  ];

  (* API *)
  SCB`Init[opts:OptionsPattern[]] := Module[{missing,PythonUtilitiesPath,defaultStatePath},
    missing = Select[{GammaLoopPath}, OptionValue[#] === Automatic &];
    If[missing =!= {}, Message[Init::req, Row[missing, ", "]]; Return[$Failed]];
    
    InfoMessage["Initializing SymbolicaCommunityBridge..."];
    
    $initialized = False;
    
    Catch[
    
    SetEnvironment["PYTHONPATH" -> ""];
    If[ Not[OptionValue[FORMPath] === Automatic],
       SetEnvironment["PATH" -> DirectoryName[OptionValue[FORMPath]]<>":"<>GetEnvironment["PATH"][[2]]];
    ];

    If[
      Not[OptionValue[SymbolicaCommunityPath] === Null],
      If[ Not[OptionValue[SymbolicaCommunityPath] === Automatic],
      SetEnvironment["PYTHONPATH" -> OptionValue[SymbolicaCommunityPath]<>":"<>GetEnvironment["PYTHONPATH"][[2]]]
      ];
    ];
    
    If[
      Not[OptionValue[SymbolicaLicense] === Automatic],
      SetEnvironment["SYMBOLICA_LICENSE" -> OptionValue[SymbolicaLicense]]
    ];

    If[
      Not[OptionValue[SymbolicaCommunityPath] === Null],
      If[ Not[OptionValue[SymbolicaCommunityPath] === Automatic],
      SetEnvironment["PYTHONPATH" -> OptionValue[SymbolicaCommunityPath]<>":"<>GetEnvironment["PYTHONPATH"][[2]]]
      ];
    ];

    If[
      Not[OptionValue[GammaLoopPath] === Null],
      If[ Not[OptionValue[GammaLoopPath] === Automatic],
        SetEnvironment["PYTHONPATH" -> OptionValue[GammaLoopPath]<>":"<>GetEnvironment["PYTHONPATH"][[2]]]
      ];
	  If[ Not[OptionValue[DebugLevel] === Automatic],
	     SetEnvironment["GL_NO_HARD_WARNINGS" -> "on"];
	     SetEnvironment["GL_DISPLAY_FILTER" -> OptionValue[DebugLevel]];
	     SetEnvironment["GL_LOGFILE_FILTER" -> OptionValue[DebugLevel]];
	  ];
    ];

    PythonUtilitiesPath = ParentDirectory[DirectoryName[SCB`$ThisFile]];
    SetEnvironment[ "PYTHONPATH" -> PythonUtilitiesPath <> ":" <> GetEnvironment["PYTHONPATH"][[2]] ];

    If[
      Not[OptionValue[PythonInterpreter] === Automatic],
      RegisterExternalEvaluator["Python",OptionValue[PythonInterpreter]];
    ];

    If[ Length[FindExternalEvaluators["Python"]] == 0,
    Message[Cmd::error, "Could not find external Python evaluator. You can specify it with the option PythonInterpreter of SCB`Init."];
    Return[$Failed];
    ];
	
	InfoMessage["Ignore deprecation warnings appearing below, if any"];
    py = If[Not[OptionValue[PythonInterpreter] === Automatic],
      StartExternalSession[<|"System" -> "Python","ReturnType" -> "Expression", "Evaluator" -> OptionValue[PythonInterpreter]|>]
      ,
      StartExternalSession[<|"System" -> "Python", "ReturnType" -> "Expression"|>]
     ];
    
    If[
      Not[OptionValue[SymbolicaCommunityPath] === Null],
      SCB`WrappedExternalEvaluate[ py, "
from symbolica import *
from symbolica.community.idenso import *
from symbolica.community.spenso import Tensor, TensorName as N, LibraryTensor, TensorNetwork,Representation,TensorStructure,TensorIndices,Tensor,Slot,TensorLibrary, ExecutionMode
"];
    ];

    If[
      Not[OptionValue[GammaLoopPath] === Null],       
      SCB`WrappedExternalEvaluate[ py, "from gammaloop import GammaLoopAPI, LogLevel"];
      If[
       OptionValue[GammaLoopStatePath] === Automatic,
       defaultStatePath = FileNameJoin[NotebookDirectory[],"gammaloop_state"];
       If[DirectoryQ[defaultStatePath],
         DeleteDirectory[defaultStatePath, DeleteContents -> True]
       ];
       ,
       defaultStatePath = OptionValue[GammaLoopStatePath];
      ];
      SCB`WrappedExternalEvaluate[ py, "gl = GammaLoopAPI(\"" <> defaultStatePath <> "\",log_file_name=\"mathematica\")"];
    ];

    (* Print[GetEnvironment["PYTHONPATH"]]; *)
    SCB`WrappedExternalEvaluate[ py, "import SCBPython as scb"];
    
    (* ExternalEvaluate[ py, "scb.hello_world()"]; *)

    $config = <|
       "PythonInterpreter" -> OptionValue[PythonInterpreter],
       "GammaLoopPath" -> OptionValue[GammaLoopPath], 
       "SymbolicaCommunityPath" -> OptionValue[SymbolicaCommunityPath], 
       "SymbolicaLicense" -> OptionValue[SymbolicaLicense],
       "ParsingModeFullForm" -> OptionValue[ParsingModeFullForm],
       "DebugLevel" -> If[OptionValue[DebugLevel]==="Info",0,1],
       "py" -> py
    |>;
     
    $initialized = True;
    
    InfoMessage["SymbolicaCommunityBridge successfully initialized !"];
    
    Return[True];
    ,
      "PythonError"
    ];
    Return[False];
  ];
  
  SCB`SetDebugLevel[DebugLevel_] := Block[{},
    $config["DebugLevel"] = DebugLevel;
  ];
  
  SCB`InitializedQ[] := TrueQ[$initialized];
  SCB`Reset[] := ($initialized = False; $config = <||>; Null);

  (* implementations *)
  (*
  SCB`SExpand[expr_] := Block[{ inputExpr, pythoncmd, res },
     
     inputExpr = SCB`ExprToString[expr, BackTickReplace -> True];
     
     pythoncmd = StringTemplate["
E(r'''`inputExpr`'''.replace('MATHEMATICABACKTICK','`'), ParseMode.Mathematica).expand().to_mathematica()
     "][<|"inputExpr"\[Rule]inputExpr|>];
     
     (* Print[pythoncmd]; *)
     
     res = WrappedExternalEvaluate[ $config["py"] , pythoncmd ];
     
     SCB`StringToExpr[res]
  ];
  *)
  
  SCB`ExpressionOrString[e_]:=Block[{},
     If[StringQ[e],
       e,
       StringTemplate["scb.to_symbolica(r'''`inputExpr`''', debug=`debug`)"][<|"inputExpr" -> SCB`ExprToString[e,BackTickReplace->True], "debug" -> If[$config["DebugLevel"]>0,"True","False"]|>]
     ]
  ];
  
  SCB`Run[expr_, cmd_?StringQ, OptionsPattern[{args -> {}, opts -> <||>}]] := Block[{pythoncmd, inputExpr, argsStr, optsStr, argsAndOptsString, stringCmd},
    
     inputExpr = SCB`ExprToString[expr, BackTickReplace -> True];
     stringCmd = ToString[cmd];
     argsStr = StringRiffle[Table[SCB`ExpressionOrString[a], {a, OptionValue[args]}],", "];
     optsStr = StringRiffle[Table[ToString[k]<>"="<>SCB`ExpressionOrString[OptionValue[opts][k]], {k, Keys[OptionValue[opts]]}],", "];
     argsAndOptsString = "";
     If[Length[OptionValue[args]]>0,
       argsAndOptsString = argsAndOptsString <> argsStr;
     ];
     If[Length[OptionValue[opts]]>0,
       argsAndOptsString = argsAndOptsString <> If[Length[OptionValue[args]]>0, ", ", ""] <>optsStr;
     ];
     pythoncmd = StringTemplate["
scb.to_mathematica_form(scb.to_symbolica(r'''`inputExpr`''', debug=`debug`).`stringcmd`(`argsandopts`).to_mathematica())
     "][<|
        "inputExpr" -> inputExpr,
        "debug" -> If[$config["DebugLevel"]>0,"True","False"],
        "stringcmd" -> stringCmd,
        "argsandopts" -> argsAndOptsString
       |>];
       
     WrappedExternalEvaluate[ $config["py"] , pythoncmd ]
  ];

  SCB`RunCmds[exprs_?AssociationQ, cmds_?StringQ] := Block[{definitions, definitionsList, pythoncmd, inputExpr},
    definitionsList = Table[ToString[k]<>" = "<>SCB`ExpressionOrString[exprs[k]],{k, Keys[exprs]}];
    definitions = StringRiffle[definitionsList, "\n"];
     pythoncmd = StringTemplate["
`inputExprs`

# Initialize conventional name for result to a default value 
res = None;

`cmds`

scb.to_mathematica_form(res.to_mathematica()) if res is not None else 'Null'

"][<|
        "inputExprs" -> definitions,
        "cmds" -> cmds,
        "debug" -> If[$config["DebugLevel"]>0,"True","False"]
       |>];
       
     WrappedExternalEvaluate[ $config["py"] , pythoncmd ]
  ];

  SCB`SExpand[expr_] := Block[{ inputExpr, pythoncmd, res },
     
     inputExpr = SCB`ExprToString[expr, BackTickReplace -> True];
     
     pythoncmd = StringTemplate["
e = scb.to_symbolica(r'''`inputExpr`''', debug=`debug`)
e = e.expand()
scb.to_mathematica_form(e.to_mathematica())
     "][<|"inputExpr"->inputExpr, "debug" -> If[$config["DebugLevel"]>0,"True","False"]|>];
     
     (* Print[pythoncmd]; *)
     
     res = WrappedExternalEvaluate[ $config["py"] , pythoncmd ];
     
     SCB`StringToExpr[res]
  ];

  SCB`SimplifyColor[expr_] := Block[{ inputExpr, pythoncmd, res },
     
     inputExpr = SCB`ExprToString[expr, BackTickReplace -> True];
     
     pythoncmd = StringTemplate["
e = scb.to_symbolica(r'''`inputExpr`''', debug=`debug`)
e = simplify_color(e)
scb.to_mathematica_form(e.to_mathematica())
     "][<|"inputExpr"->inputExpr, "debug" -> If[$config["DebugLevel"]>0,"True","False"]|>];
     
     (* Print[pythoncmd]; *)
     
     res = WrappedExternalEvaluate[ $config["py"] , pythoncmd ];
     
     SCB`StringToExpr[res]
  ];
  
  SCB`SimplifyGamma[expr_] := Block[{ inputExpr, pythoncmd, res },
     
     inputExpr = SCB`ExprToString[expr, BackTickReplace -> True];
     
     pythoncmd = StringTemplate["
e = scb.to_symbolica(r'''`inputExpr`''', debug=`debug`)
e = simplify_gamma(e)
scb.to_mathematica_form(e.to_mathematica())
     "][<|"inputExpr"->inputExpr, "debug" -> If[$config["DebugLevel"]>0,"True","False"]|>];
     
     (* Print[pythoncmd]; *)
     
     res = WrappedExternalEvaluate[ $config["py"] , pythoncmd ];
     
     SCB`StringToExpr[res]
  ];
  
  SCB`RunPython[pythonCode_] := Block[{},
    WrappedExternalEvaluate[ $config["py"] , pythonCode]
  ];

  SCB`ResetGammaLoop[OptionsPattern[{GammaLoopStatePath -> Automatic}]] := Module[
     {}, 
     If[
       OptionValue[GammaLoopStatePath] === Automatic,
       defaultStatePath = FileNameJoin[NotebookDirectory[],"gammaloop_state"];
       If[DirectoryQ[defaultStatePath],
         DeleteDirectory[defaultStatePath, DeleteContents -> True]
       ];
       ,
       defaultStatePath = OptionValue[GammaLoopStatePath];
      ];
      SCB`WrappedExternalEvaluate[ py, "gl = GammaLoopAPI(\"" <> defaultStatePath <> "\",log_file_name=\"mathematica\")"];
  ];
      
  SCB`GammaLoop[cmd_, opts:OptionsPattern[]] := Module[
     {lvl = OptionValue[DebugLevel]}, 

     pythoncmd = StringTemplate["
gl.run(\"`cmd`\")
     "][<| "cmd" -> cmd |>];
     
     res = WrappedExternalEvaluate[ $config["py"] , pythoncmd ];
     
     res

  ];
  
  SCB`GammaLoopImportGraphs[graphs_?StringQ, processName_?StringQ, OptionsPattern[{Overwrite -> True}]] := Module[{},

     pythoncmd = StringTemplate["
gl.import_graphs(\"\"\"`graphs`\"\"\",process_name=\"`process_name`\",format='string',overwrite=`do_overwrite`)
     "][<| "graphs" -> graphs, "process_name" -> processName, "do_overwrite" -> If[OptionValue[Overwrite], "True", "False"] |>];

     res = WrappedExternalEvaluate[ $config["py"] , pythoncmd ];
  
  ];
  
];

(* one entry point to (re)define everything; manages protect/unprotect *)
iLoad[mode_] := Module[{},
  iUnprotectPublic[];
  iClearImpl[];
  iDefineImpl[];
  iProtectPublic[];
  If[mode === "init",
    Message[SCB`Welcome::info]
  ];
];

(* first-time attach *)
iLoad[If[TrueQ[Global`SCB$Reloading] && NameQ["SCB`Private`iLoad"],"reload","init"]];

End[];  (* `Private` *)

(* public hot-reload: re-Get this file; the top trap runs iLoad[] and aborts further eval *)
SCB`LoadSCB[] := Module[{},
  Global`SCB$Reloading = True;
  Catch[
    Block[{$ContextPath = {"System`"}}, Get[SCB`$ThisFile]],
    "SCB/Reload"
  ];
  Global`SCB$Reloading = False;
  True
];

EndPackage[];

(* force qualification: keep SCB` off $ContextPath *)
If[MemberQ[$ContextPath, "SCB`"], $ContextPath = DeleteCases[$ContextPath, "SCB`"]];
