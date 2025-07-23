{
  description = "Naive Nix packaging for KGD";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

    flake-utils.url = "github:numtide/flake-utils";
  };

  # KGD is not a first-class R package, merely a collection of R source files
  # to be source'd into the caller's environment.
  #
  # The Nix package enables these R sources to be used directly from the Nix store.

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
          };

          KGD-rPackages = with pkgs.rPackages;
            [
              RcppArmadillo
              data_table
              R_utils
              plotly
              heatmaply
              parallelDist
              MethComp
              MASS
            ];

          KDG-src = pkgs.stdenv.mkDerivation {
            pname = "KDG-src";
            version = "1.3.1";

            src = ./.;

            nativeBuildInputs = [ pkgs.dos2unix ];

            installPhase = ''
              mkdir $out
              runHook preInstall
              cp GBS-Chip-Gmatrix.R GBSPedAssign.R GBS-PopGen.R GBSRun.R GBS-Rcpp-functions.cpp $out
              chmod 755 $out/GBSRun.R
              dos2unix $out/*
              runHook postInstall
            '';

            postFixup = ''
              substituteInPlace $out/GBSRun.R --replace '<source directory>' $out
            '';
          };
        in
        {
          packages = {
            src = KDG-src;
            rPackages = KGD-rPackages;
          };
        });
}
