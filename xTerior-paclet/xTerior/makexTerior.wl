(* ::Package:: *)

PacletUninstall["xTerior"];
xTerior = CreatePacletArchive[DirectoryName[$InputFileName]];
PacletInstall[xTerior];
