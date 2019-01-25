.onAttach <-
function(libname, pkgname) {
    packageStartupMessage("\nBENMMI version: ", utils::packageVersion("BENMMI"))
    packageStartupMessage(
        "Copyright 2015-", 
        format(Sys.Date(), "%Y"), 
        " by Rijkswaterstaat, the Netherlands (RWS)."
    )
    packageStartupMessage("Type  citation(\"BENMMI\")  on how to cite BENMMI in publications.")
    packageStartupMessage("For the tutorial, type: vignette(\"BENMMI_User_Manual\")")
    packageStartupMessage("For Frequently Asked Questions, type: vignette(\"FAQ\")")

}