
Module[{root = FileNameJoin[{DirectoryName @ $InputFileName, ".."}],
        file},
  file = FileNameJoin[{root, "SCB", "SCB.wl"}];

  Quiet[
    Block[{$ContextPath = {"System`"}},      (* keep Global` off the path during load *)
      Get[file]
    ],
    General::shdw                           (* suppress any residual shadowing from deps *)
  ];
]

