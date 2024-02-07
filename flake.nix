{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
    devenv.url = "github:cachix/devenv";
  };

  nixConfig = {
    extra-trusted-public-keys = "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };

  outputs = { self, nixpkgs, devenv, systems, ... } @ inputs:
    let
      forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in {
      packages = forEachSystem (system: 
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in rec {
        devenv-up = self.devShells.${system}.default.config.procfileScript;
      });

      devShells = forEachSystem (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in rec {
          default = devenv.lib.mkShell {
            inherit inputs pkgs;
            modules = [{
                packages = with pkgs; with python3Packages; [
                  gnumake
                  numpy
                  matplotlib
                  packaging
                  pillow
                  pyparsing
                  cycler
                  dateutil
                  kiwisolver
                  six
                ];
                
                languages.cplusplus.enable = true;
                languages.python.enable = true;
            }];
          };
        });
    };
}

