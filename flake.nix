{
  description = "Naive Nix packaging for KGD";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";

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

          # we need a recent Rcpp with this fix:
          # https://github.com/RcppCore/Rcpp/pull/1346
          recent-Rcpp = pkgs.rPackages.buildRPackage {
            name = "Rcpp";

            version = "1.0.13.6";

            src = pkgs.fetchFromGitHub {
              owner = "RcppCore";
              repo = "Rcpp";
              rev = "83e640b55aeaeba8a746d0da6152cabe8af41154";
              hash = "sha256-x0s9BRIAxxZmNFXRYwtpmT5asD4+mhIUTqwYsaPOgi4=";
            };
          };

          # TODO: it should be possible to simply override the Rcpp dependency for RcppArmadillo in nixpkgs,
          # but I didn't work out how to do that
          RcppArmadillo-with-recent-Rcpp = pkgs.rPackages.buildRPackage {
            name = "RcppArmadillo";

            version = "0.12.8.1.0";

            src = pkgs.fetchFromGitHub {
              owner = "RcppCore";
              repo = "RcppArmadillo";
              rev = "8abe7be9fc4dd7c1d2b02ed200707232d6fd1f09"; # 0.12.8.1.0
              hash = "sha256-+Li4ln/4ZyBY+I8S8X4uSmFaG1D3q5UJOJJB5pLRubo=";
            };

            propagatedBuildInputs = [ recent-Rcpp ];
          };


          KDG-R = pkgs.rWrapper.override
            {
              packages = with pkgs.rPackages; [
                RcppArmadillo-with-recent-Rcpp
                data_table
                R_utils
                plotly
                heatmaply
                parallelDist
                MethComp
                MASS
              ];
            };

          KDG-src = pkgs.stdenv.mkDerivation {
            pname = "KDG-src";
            version = "1.2.2";

            src = ./.;

            buildInputs = [ KDG-R ];

            propagatedBuildInputs = [ KDG-R ];

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
          };
        });
}
