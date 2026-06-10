{
  description = "rabe pinned build/test toolchain for reproducible golden records";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "aarch64-darwin" "x86_64-darwin" ];
      forAll = f: nixpkgs.lib.genAttrs systems (s: f nixpkgs.legacyPackages.${s});
      # Toolchain pinned to nixpkgs so gfortran, the C/Fortran runtime, libm,
      # cmake and netcdf are byte-identical on any machine using this flake.
      # libneo is left to the in-tree FetchContent pin (cmake/libneo.cmake).
      toolchain = pkgs: [
        pkgs.gfortran
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
    in {
      devShells = forAll (pkgs: {
        default = pkgs.mkShell { packages = toolchain pkgs ++ [ (pyenv pkgs) ]; };
        # Toolchain only: substitutable binaries, no local build step, so it works
        # under proot / no-user-namespace sandboxes (HPC login nodes, Condor
        # execute nodes). Run the python compare with the host interpreter.
        build = pkgs.mkShell { packages = toolchain pkgs; };
      });
    };
}
