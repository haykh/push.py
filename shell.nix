{
  pkgs ? import <nixpkgs> { },
  py ? "312",
}:

pkgs.mkShell {
  name = "pushpydev";
  nativeBuildInputs = with pkgs; [
    pkgs."python${py}"
    pkgs."python${py}Packages".pip
    black
    pyright
    taplo
    vscode-langservers-extracted
    zlib
  ];
  LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath [
    pkgs.stdenv.cc.cc
    pkgs.zlib
  ];

  shellHook = ''
    if [ ! -d ".venv" ]; then
      ./deploy.sh
    else
      source .venv/bin/activate
    fi
    echo "push.py nix-shell activated: $(which python)"
  '';
}
