.onLoad <- function(lib, pkgname) {
  if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("FunciSNP.data")
    }
}

