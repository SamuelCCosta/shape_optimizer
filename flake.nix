{
description = "Neoseminario + python time testing";

inputs = {
    nixpkgs.url = "github:Nixos/nixpkgs/nixos-25.05";
};

outputs = { self, nixpkgs }: {
    devShells.x86_64-linux.default =
        let
            system = "x86_64-linux";
            pkgs = nixpkgs.legacyPackages.${system};

            maniFEM-cflags-list = [
                "-O3"
                "-std=c++23"
                "-c"
                "-Wshadow"
                "-Wall"
                "-fPIC"
                "-I ."
                "-I ${pkgs.eigen}/include/eigen3"
                "-DOMIT_OBSOLETE_CODE"

                # You may comment or uncomment lines below
                # "-DMANIFEM_NO_FEM"
                # "-DMANIFEM_NO_FRONTAL"
                # "-DMANIFEM_NO_QUOTIENT"
                # "-DMANIFEM_COLLECT_CM"
                # "-DNDEBUG"
            ];

            maniFEM-cflags = pkgs.lib.concatStringsSep " " maniFEM-cflags-list;

            # Derivation for the maniFEM library
            maniFEM = pkgs.stdenv.mkDerivation rec {
                pname = "maniFEM";
                version = "25.02";

                src = pkgs.fetchgit {
                    url = "https://codeberg.org/cristian.barbarosie/maniFEM.git";
                    rev = "ca6d32bcdf7c35ae8d03a1680d40ee4ff9e643b3";
                    sha256 = "sha256-FVNobCSV01llWUbfVp+86mVDQK9qbat+4rxIjMnc0Do=";
                };

                buildInputs = [ pkgs.eigen ];

                buildPhase = ''
                    runHook preBuild
                    make static-lib "CFLAGS=${maniFEM-cflags}"
                    runHook postBuild
                '';

                installPhase = ''
                    runHook preInstall
                    install -d $out/lib $out/include
                    install -m644 libmaniFEM.a $out/lib
                    install -m644 maniFEM.h $out/include
                    install -d $out/include/src
                    install -m644 src/*.h $out/include/src
                    runHook postInstall
                '';
            };

            python-with-packages = pkgs.python3.withPackages (ps: [
                ps.pandas
                ps.pybind11
                ps.pybind11-stubgen
                ps.matplotlib
            ]);

        in
        pkgs.mkShell {
            packages = [
                pkgs.gcc
                pkgs.gnumake
                pkgs.eigen
                maniFEM
                python-with-packages
            ];
          
            CPLUS_INCLUDE_PATH = "${maniFEM}/include:${pkgs.eigen}/include/eigen3";
            LIBRARY_PATH = "${maniFEM}/lib";

            shellHook = ''
                export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$(python3 -c "import pybind11; print(pybind11.get_include())")
            '';

            NIX_CFLAGS_COMPILE = "-I${maniFEM}/include";
            NIX_LDFLAGS = "-L${maniFEM}/lib";
        };
};
}