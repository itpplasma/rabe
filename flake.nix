{
  description = "rabe pinned build/test toolchain for reproducible golden records";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/d6133526472eb11a863eb6e679b104086ef291bf";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" ];
      forAll = f: nixpkgs.lib.genAttrs systems (s: f nixpkgs.legacyPackages.${s});
      toolchain = pkgs: [
        pkgs.gcc14
        pkgs.gfortran14
        pkgs.cmake
        pkgs.gnumake
        pkgs.ninja
        pkgs.pkg-config
        pkgs.git
        pkgs.netcdf
        pkgs.netcdffortran
      ];
      # Python for the golden-record compare; only needed where the test runs.
      pyenv = pkgs: pkgs.python3.withPackages (ps: [ ps.xarray ps.numpy ps.netcdf4 ]);
      mkShell = pkgs: packages:
        (pkgs.mkShell.override { stdenv = pkgs.gcc14Stdenv; }) {
          packages = packages;
          CC = "${pkgs.gcc14}/bin/gcc";
          CXX = "${pkgs.gcc14}/bin/g++";
          FC = "${pkgs.gfortran14}/bin/gfortran";
        };
    in {
      devShells = forAll (pkgs: {
        default = mkShell pkgs (toolchain pkgs ++ [ (pyenv pkgs) ]);
      });
    };
}
