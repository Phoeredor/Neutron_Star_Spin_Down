#!/usr/bin/env wolframscript

(* Physical constants *)

ClearAll["Global`*"];

secondsPerYear = 365.25*24.*3600.;
spinDownConstant = 2.44*10^-40; (* s G^-2 *)
weakMagneticField = 1.0*10^13; (* G *)
strongMagneticField = 5.0*10^14; (* G *)
decayTimescaleYears = 1.0*10^4; (* yr *)
decayTimescaleSeconds = decayTimescaleYears*secondsPerYear;
maximumTimeYears = 1.0*10^6; (* yr *)
initialPeriods = <|"1 ms" -> 1.0*10^-3, "3 ms" -> 3.0*10^-3|>; (* s *)

stellarRadius = 1.0*10^6; (* cm *)
momentOfInertia = 1.0*10^45; (* g cm^2 *)
solarMass = 2.0*10^33; (* g *)
stellarMass = 1.25*solarMass; (* g *)
gravitationalConstant = 6.67259*10^-8; (* cm^3 g^-1 s^-2 *)
periodDerivativeCutoff = 1.0*10^-14; (* s/s *)

(* Model functions *)

constantFieldPeriod[timeSeconds_, initialPeriod_, magneticField_] :=
  Sqrt[initialPeriod^2 + 2*spinDownConstant*magneticField^2*timeSeconds];

constantFieldPeriodDerivative[timeSeconds_, initialPeriod_, magneticField_] :=
  spinDownConstant*magneticField^2/
    constantFieldPeriod[timeSeconds, initialPeriod, magneticField];

decayingFieldPeriod[
    timeSeconds_, initialPeriod_, initialMagneticField_, timescaleSeconds_] :=
  Sqrt[
    initialPeriod^2 +
      spinDownConstant*initialMagneticField^2*timescaleSeconds*
        (1 - Exp[-2*timeSeconds/timescaleSeconds])
  ];

decayingFieldPeriodDerivative[
    timeSeconds_, initialPeriod_, initialMagneticField_, timescaleSeconds_] :=
  spinDownConstant*initialMagneticField^2*Exp[-2*timeSeconds/timescaleSeconds]/
    decayingFieldPeriod[
      timeSeconds, initialPeriod, initialMagneticField, timescaleSeconds];

angularVelocity[periodSeconds_] := 2*Pi/periodSeconds;

relativePercentageVariation[periodValues_, initialPeriod_] :=
  100*(periodValues - initialPeriod)/initialPeriod;

normalizedPeriodDerivative[derivativeValues_] :=
  derivativeValues/First[derivativeValues];

asymptoticDecayingFieldPeriod[
    initialPeriod_, initialMagneticField_, timescaleSeconds_] :=
  Sqrt[
    initialPeriod^2 +
      spinDownConstant*initialMagneticField^2*timescaleSeconds
  ];

(* Derived quantities *)

fieldPrescriptions = <|
  "Constant weak field" -> <|
    "Model" -> "Constant", "Field" -> weakMagneticField|>,
  "Constant strong field" -> <|
    "Model" -> "Constant", "Field" -> strongMagneticField|>,
  "Exponentially decaying strong field" -> <|
    "Model" -> "Decaying", "Field" -> strongMagneticField|>
|>;

computeTrack[initialPeriod_, specification_Association] := Module[
  {periodValues, derivativeValues, asymptoticPeriod},
  If[specification["Model"] === "Constant",
    periodValues = constantFieldPeriod[
      timeSeconds, initialPeriod, specification["Field"]];
    derivativeValues = constantFieldPeriodDerivative[
      timeSeconds, initialPeriod, specification["Field"]];
    asymptoticPeriod = Missing["NotApplicable"],
    periodValues = decayingFieldPeriod[
      timeSeconds, initialPeriod, specification["Field"],
      decayTimescaleSeconds];
    derivativeValues = decayingFieldPeriodDerivative[
      timeSeconds, initialPeriod, specification["Field"],
      decayTimescaleSeconds];
    asymptoticPeriod = asymptoticDecayingFieldPeriod[
      initialPeriod, specification["Field"], decayTimescaleSeconds]
  ];
  <|
    "Period" -> periodValues,
    "PeriodDerivative" -> derivativeValues,
    "AngularVelocity" -> angularVelocity[periodValues],
    "RelativePercentageVariation" ->
      relativePercentageVariation[periodValues, initialPeriod],
    "NormalizedDerivative" -> normalizedPeriodDerivative[derivativeValues],
    "AsymptoticPeriod" -> asymptoticPeriod
  |>
];

computeTracks[initialPeriod_] := Association[
  KeyValueMap[
    Function[{fieldName, specification},
      fieldName -> computeTrack[initialPeriod, specification]],
    fieldPrescriptions
  ]
];

(* Scientific consistency checks *)

failCheck[message_String] := (
  Print["Scientific consistency error: " <> message];
  Exit[1]
);

ensureCheck[condition_, message_String] :=
  If[!TrueQ[condition], failCheck[message]];

finiteRealNumberQ[value_] :=
  NumericQ[value] &&
    FreeQ[value, Indeterminate | ComplexInfinity | DirectedInfinity | _Complex];

finiteRealVectorQ[values_] := VectorQ[values, finiteRealNumberQ];

validateTrack[
    label_String, initialPeriod_, track_Association, decayingModel_] := Module[
  {periodValues, derivativeValues, generatedArrays, monotonicTolerance,
   asymptoticPeriod},
  periodValues = track["Period"];
  derivativeValues = track["PeriodDerivative"];
  generatedArrays = {
    periodValues, derivativeValues, track["AngularVelocity"],
    track["RelativePercentageVariation"], track["NormalizedDerivative"]
  };
  ensureCheck[
    And @@ (finiteRealVectorQ /@ generatedArrays),
    label <> " contains a non-finite or non-real value."];
  ensureCheck[And @@ Thread[periodValues > 0],
    label <> " contains a non-positive period."];
  ensureCheck[And @@ Thread[derivativeValues >= 0],
    label <> " contains a negative period derivative."];
  monotonicTolerance = 64*$MachineEpsilon*Max[periodValues];
  ensureCheck[Min[Differences[periodValues]] >= -monotonicTolerance,
    label <> " has a period that decreases with time."];
  ensureCheck[Abs[First[periodValues] - initialPeriod] <= 10^-15,
    label <> " does not begin at the requested initial period."];
  If[TrueQ[decayingModel],
    asymptoticPeriod = track["AsymptoticPeriod"];
    ensureCheck[
      Max[periodValues] <= asymptoticPeriod*(1 + 10^-12),
      label <> " exceeds the analytical decaying-field asymptote."];
    ensureCheck[
      Abs[Last[periodValues] - asymptoticPeriod] <=
        10^-10*Max[1, asymptoticPeriod],
      label <> " disagrees with the analytical asymptote at late time."]
  ]
];

(* Centrifugal-stability estimate *)

breakupAngularVelocity[gravity_, mass_, radius_] := Sqrt[gravity*mass/radius^3];
breakupPeriod[gravity_, mass_, radius_] :=
  2*Pi/breakupAngularVelocity[gravity, mass, radius];

breakupAngularVelocityValue = breakupAngularVelocity[
  gravitationalConstant, stellarMass, stellarRadius];
breakupPeriodValue = breakupPeriod[
  gravitationalConstant, stellarMass, stellarRadius];

ensureCheck[
  And @@ Thread[Values[initialPeriods] > breakupPeriodValue],
  "An initial period is below the centrifugal breakup estimate."];

(* Time grid *)

positiveTimeYears = 10.^Subdivide[-6., Log10[maximumTimeYears], 3000];
timeYears = Join[{0.}, positiveTimeYears];
timeSeconds = timeYears*secondsPerYear;

ensureCheck[finiteRealVectorQ[timeYears], "The time grid is not finite."];
ensureCheck[Min[Differences[timeYears]] > 0,
  "The time grid is not strictly increasing."];

allResults = Association[
  KeyValueMap[
    Function[{initialLabel, initialPeriod},
      initialLabel -> computeTracks[initialPeriod]],
    initialPeriods
  ]
];

KeyValueMap[
  Function[{initialLabel, tracks},
    KeyValueMap[
      Function[{fieldName, track},
        validateTrack[
          initialLabel <> ", " <> fieldName,
          initialPeriods[initialLabel],
          track,
          fieldPrescriptions[fieldName]["Model"] === "Decaying"
        ]
      ],
      tracks
    ]
  ],
  allResults
];

(* Figure construction *)

trackNames = Keys[fieldPrescriptions];
fieldPlotStyles = <|
  "Constant weak field" -> Directive[RGBColor[0.19, 0.41, 0.56], Thick],
  "Constant strong field" ->
    Directive[RGBColor[0.71, 0.21, 0.47], Dashed, Thick],
  "Exponentially decaying strong field" ->
    Directive[RGBColor[0.21, 0.72, 0.48], DotDashed, Thick]
|>;

(* ArcSinh is a smooth symmetric-logarithmic coordinate with a valid origin. *)
timeLinearScaleYears = 10^-3;
timeCoordinate[value_] := ArcSinh[value/timeLinearScaleYears];
timeCoordinates = timeCoordinate /@ timeYears;
timeTickValues = {0., 10^-3, 10^-2, 10^-1, 1., 10., 10^2, 10^3,
  10^4, 10^5, 10^6};
timeTicks = ({timeCoordinate[#], If[# == 0., "0", ScientificForm[#, 1]]} & /@
  timeTickValues);

makeTimeSeriesPlot[
    tracks_Association, quantity_String, verticalLabel_String,
    logarithmicVertical_] := ListLinePlot[
  Table[
    Transpose[{timeCoordinates, tracks[fieldName][quantity]}],
    {fieldName, trackNames}
  ],
  Joined -> True,
  PlotStyle -> Values[fieldPlotStyles],
  PlotLegends -> Placed[LineLegend[Values[fieldPlotStyles], trackNames], Below],
  ScalingFunctions -> If[TrueQ[logarithmicVertical], {None, "Log10"}, None],
  Frame -> True,
  FrameTicks -> {{Automatic, Automatic}, {timeTicks, None}},
  FrameLabel -> {{verticalLabel, None}, {"Time [yr]", None}},
  GridLines -> Automatic,
  PlotRange -> All,
  ImageSize -> 520
];

makeSpinEvolutionFigure[initialLabel_String, tracks_Association] := Labeled[
  Grid[{
    {
      makeTimeSeriesPlot[tracks, "RelativePercentageVariation",
        "Relative period variation [%]", False],
      makeTimeSeriesPlot[tracks, "Period", "Period [s]", True]
    },
    {
      makeTimeSeriesPlot[tracks, "PeriodDerivative",
        "Period derivative [s/s]", True],
      makeTimeSeriesPlot[tracks, "NormalizedDerivative",
        "Normalized derivative Pdot(t)/Pdot(0)", True]
    }
  }, Spacings -> {1.2, 1.2}],
  Style["Neutron-star spin-down evolution (P0 = " <> initialLabel <> ")", 18],
  Top
];

initialPeriodColors = <|
  "1 ms" -> RGBColor[0.23, 0.32, 0.55],
  "3 ms" -> RGBColor[0.91, 0.44, 0.32]
|>;

(* A smaller initial period produces a larger normalized percentage change;
   the absolute constant-field periods progressively converge at late time. *)
initialPeriodComparisonFigure = ListLinePlot[
  Table[
    Transpose[{
      timeCoordinates,
      allResults[initialLabel]["Constant strong field"][
        "RelativePercentageVariation"]
    }],
    {initialLabel, Keys[initialPeriods]}
  ],
  Joined -> True,
  PlotStyle -> (Directive[#, Thick] & /@ Values[initialPeriodColors]),
  PlotLegends -> Placed[
    LineLegend[Values[initialPeriodColors],
      ("P0 = " <> # & /@ Keys[initialPeriods])], Below],
  Frame -> True,
  FrameTicks -> {{Automatic, Automatic}, {timeTicks, None}},
  FrameLabel -> {{"Relative period variation [%]", None}, {"Time [yr]", None}},
  PlotLabel -> "Initial-period comparison for constant B = 5 x 10^14 G",
  GridLines -> Automatic,
  PlotRange -> All,
  ImageSize -> 1000
];

fieldDashing = <|
  "Constant weak field" -> Dashing[{}],
  "Constant strong field" -> Dashing[{0.025, 0.015}],
  "Exponentially decaying strong field" -> Dashing[{0.025, 0.01, 0.005, 0.01}]
|>;

periodDerivativeSeries = Flatten[
  Table[
    Select[
      Transpose[{
        allResults[initialLabel][fieldName]["Period"],
        allResults[initialLabel][fieldName]["PeriodDerivative"]
      }],
      Last[#] >= periodDerivativeCutoff &
    ],
    {initialLabel, Keys[initialPeriods]}, {fieldName, trackNames}
  ],
  1
];

ensureCheck[And @@ Thread[Length /@ periodDerivativeSeries > 0],
  "A P-Pdot trajectory has no samples above the adopted cutoff."];

periodDerivativeStyles = Flatten[
  Table[
    Directive[
      initialPeriodColors[initialLabel], fieldDashing[fieldName],
      AbsoluteThickness[2]
    ],
    {initialLabel, Keys[initialPeriods]}, {fieldName, trackNames}
  ],
  1
];

periodDerivativeLabels = Flatten[
  Table[
    initialLabel <> " - " <> ToLowerCase[fieldName],
    {initialLabel, Keys[initialPeriods]}, {fieldName, trackNames}
  ],
  1
];

displayedPeriods = Flatten[periodDerivativeSeries[[All, All, 1]]];
cutoffReferenceLine = {
  {Min[displayedPeriods], periodDerivativeCutoff},
  {Max[displayedPeriods], periodDerivativeCutoff}
};

periodDerivativeDiagram = ListLogLogPlot[
  Append[periodDerivativeSeries, cutoffReferenceLine],
  Joined -> True,
  PlotStyle -> Append[
    periodDerivativeStyles,
    Directive[GrayLevel[0.35], Dotted, AbsoluteThickness[1.5]]
  ],
  PlotLegends -> Placed[
    LineLegend[
      Append[periodDerivativeStyles,
        Directive[GrayLevel[0.35], Dotted, AbsoluteThickness[1.5]]],
      Append[periodDerivativeLabels,
        "Adopted spin-down cutoff: Pdot = 10^-14 s/s"]
    ],
    Below
  ],
  Frame -> True,
  FrameLabel -> {{"Period derivative Pdot [s/s]", None}, {"Period P [s]", None}},
  PlotLabel -> "Neutron-star evolution in the P-Pdot plane",
  GridLines -> Automatic,
  PlotRange -> All,
  ImageSize -> 1100
];

spinEvolutionOneMillisecondFigure = makeSpinEvolutionFigure[
  "1 ms", allResults["1 ms"]];
spinEvolutionThreeMillisecondFigure = makeSpinEvolutionFigure[
  "3 ms", allResults["3 ms"]];

(* Export *)

sourceDirectory = DirectoryName[$InputFileName];
repositoryRoot = DirectoryName[sourceDirectory];
outputDirectory = FileNameJoin[{repositoryRoot, "plots"}];
If[!DirectoryQ[outputDirectory], CreateDirectory[outputDirectory]];

figures = <|
  "spin_evolution_1ms.png" -> spinEvolutionOneMillisecondFigure,
  "spin_evolution_3ms.png" -> spinEvolutionThreeMillisecondFigure,
  "initial_period_comparison.png" -> initialPeriodComparisonFigure,
  "p_pdot_diagram.png" -> periodDerivativeDiagram
|>;

exportedPaths = KeyValueMap[
  Function[{fileName, figure},
    With[{path = FileNameJoin[{outputDirectory, fileName}]},
      Export[path, figure, "PNG", ImageResolution -> 300];
      path
    ]
  ],
  figures
];

(* Numerical summary *)

Print["Neutron-star spin-down summary"];
Print["  seconds per year: ", NumberForm[secondsPerYear, {12, 0}], " s"];
Print["  K: ", ScientificForm[spinDownConstant, 4], " s G^-2"];
Print["  weak magnetic field: ", ScientificForm[weakMagneticField, 4], " G"];
Print["  strong magnetic field: ", ScientificForm[strongMagneticField, 4], " G"];
Print["  decay timescale: ", ScientificForm[decayTimescaleYears, 4], " yr"];
Print["  maximum evolution time: ", ScientificForm[maximumTimeYears, 4], " yr"];
Print["  stellar radius: ", ScientificForm[stellarRadius, 4], " cm"];
Print["  moment of inertia: ", ScientificForm[momentOfInertia, 4], " g cm^2"];
Print["  stellar mass: ", ScientificForm[stellarMass, 4], " g"];
Print["  P-Pdot display cutoff: ", ScientificForm[periodDerivativeCutoff, 4], " s/s"];

KeyValueMap[
  Function[{initialLabel, initialPeriod},
    Print["  initial period ", initialLabel, ": ",
      ScientificForm[initialPeriod, 6], " s"];
    KeyValueMap[
      Function[{fieldName, track},
        Print["    ", fieldName, ", P(1 Myr): ",
          ScientificForm[Last[track["Period"]], 8], " s"]
      ],
      allResults[initialLabel]
    ];
    Print["    decaying-field asymptotic period: ",
      ScientificForm[
        allResults[initialLabel]["Exponentially decaying strong field"][
          "AsymptoticPeriod"], 8], " s"]
  ],
  initialPeriods
];

Print["  breakup angular velocity: ",
  ScientificForm[breakupAngularVelocityValue, 8], " rad/s"];
Print["  breakup period: ", ScientificForm[breakupPeriodValue, 8], " s"];
Print["Generated figures:"];
Scan[Print["  ", #] &, exportedPaths];
