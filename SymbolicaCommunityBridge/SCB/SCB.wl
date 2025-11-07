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

* Symbolica community:
    
    'python3 -m pip install symbolica'
    
    or
    
    'git clone git@github.com:benruijl/symbolica-community.git && cd symbolica-community'
    'cargo run --features \"python_stubgen\"'
    'maturin build --release'

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
SCB`$PublicSymbols = {Init, Reset, InitializedQ, SimplifyColor, SimplifyGamma, GammaLoop, RunPython, ResetGammaLoop};
SCB`$GatedSymbols  = {SimplifyColor, SimplifyGamma, GammaLoop, RunPython};

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


  Options[Init]      = {PythonInterpreter -> Automatic, GammaLoopPath -> Automatic, GammaLoopStatePath -> Automatic, SymbolicaCommunityPath -> Automatic, FORMPath -> Automatic, SymbolicaLicense -> Automatic, DebugLevel -> Automatic};
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
    res = ExternalEvaluate[py,cmd];
    If[Head[res]===Failure,
       Print[res];
       Throw[Null, "PythonError"];
    ];
    res
  ];
  
  SCB`ExprToString[expr_]:=Block[{},
     If[StringQ[expr],
       expr,
       ToString[InputForm[expr]]
     ]
  ];

  SCB`StringToExpr[expr_]:=Block[{},
     If[StringQ[expr],
       If[And[StringLength[expr]>=2, StringTake[expr,1]==="'", StringTake[expr,-1]==="'"],
         ToExpression[StringTake[expr,{2,-2}]]
         ,
         ToExpression[expr]
       ],
       expr
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
      StartExternalSession[<|"System" -> "Python","ReturnType" -> "String", "Evaluator" -> OptionValue[PythonInterpreter]|>]
      ,
      StartExternalSession[<|"System" -> "Python", "ReturnType" -> "String"|>]
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
       "DebugLevel" -> OptionValue[DebugLevel],
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

  SCB`InitializedQ[] := TrueQ[$initialized];
  SCB`Reset[] := ($initialized = False; $config = <||>; Null);

  (* implementations *)

  SCB`SimplifyColor[expr_] := 140*expr;
  
  SCB`SimplifyGamma[expr_] := Block[{ inputExpr, pythoncmd, res },
     
     inputExpr = SCB`ExprToString[expr];
     
     pythoncmd = StringTemplate["
e = E(\"`inputExpr`\")
e = simplify_gamma(e)
e.format(PrintMode.Mathematica)
     "][<|"inputExpr"->inputExpr|>];
     
     (* Print[pythoncmd];*)
     
     res = WrappedExternalEvaluate[ $config["py"] , pythoncmd ];
     
     SCB`StringToExpr[res]
  ];
  
  SCB`RunPython[pythonCode_] := Block[{},
    WrappedExternalEvaluate[ $config["py"] , pythonCode]
  ];

  SCB`ResetGammaLoop[cmd_, OptionsPattern[{GammaLoopStatePath -> Automatic}]] := Module[
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
