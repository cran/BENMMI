#' Read BENMMI Settings File
#'
#' This function reads BENMMI settings files. See the package vignette
#' for a description of its format. Create an example file by
#' calling the \code{\link[BENMMI]{BENMMIdir}}-function.
#'     
#' @param filename name of BENMMI input file (\code{character})
#'
#' @details The function performs the following tasks:
#' \itemize{
#'      \item{checks the existence of \code{filename};}
#'  	\item{reads JSON file while ignoring C-style comments;}
#'      \item{checks avaiability of required keys in the JSON-file}
#'  	\item{checks values in JSON-file}
#'  }
#'
#' @seealso \code{\link[BENMMI]{BENMMIdir}}
#' @import jsonlite
read_settings <- 
function(filename) {

    # check existence of settings file
	if (!file.exists(filename)) {
		stop("File not found", call. = FALSE)
	}

	# read settings file
	settings <- readLines(con = filename, warn = FALSE)

	# remove all C-style comments (//)
	# Note: comments are formally not part of the JSON specification!
	settings <- sub(pattern = "//.*$", replacement = "", x = settings)

	# parse JSON
	if (!validate(settings)) {
		stop(
			sprintf(
				"Errors found in %s. Check JSON-format (e.g. brackets, braces, trailing comma's)", 
				sQuote(filename)
			),  
			call. = FALSE
		)
	}
	settings <- fromJSON(settings)
    names(settings) <- tolower(names(settings))

    # define keys
    required_keys <- c("title", "user", "date", "files", "indicators",
                       "confidenceLevel", "months", "pooling", 
                       "genusToSpeciesConversion")
    optional_keys <- c("weights", "pressure", "legendText", "model")
    available_keys <- c(required_keys, optional_keys)

    # either indices or indicators are allowed, but not both
    if (!is.null(settings$indices) & !is.null(settings$indicators)) {
        stop(
            sprintf(
                fmt = "Use either key %s or key %s in JSON-file\n%s, but not both.", 
                sQuote("indices"), 
                sQuote("indicators"), 
                sQuote(basename(filename))
            ), 
            call. = FALSE
        )
    }
    
    # rename key 'indices' to 'indicators'
    if (!is.null(settings$indices)) {
        settings$indicators <- settings$indices
        settings$indices <- NULL
    }

    # check if required keys are available
    found <- tolower(required_keys) %in% names(settings)
    if (any(!found)) {
        stop(
            sprintf(
                fmt = "key %s is missing in JSON-file %s\n(see package vignette)", 
                toString(sQuote(required_keys[!found])), 
                sQuote(basename(filename))
            ), 
            call. = FALSE
        )
    }
    
    # check required files
    required_keys <- c("benthos", "taxa", "groupsToExclude", "habitats")
    names(settings$files) <- tolower(names(settings$files))
    found <- tolower(required_keys) %in% names(settings$files)
    if (any(!found)) {
        stop(
            sprintf(
                fmt = "The following keys are missing in %s: %s\n(see package vignette)", 
                sQuote(basename(filename)),
                toString(sQuote(required_keys[!found])) 
            ), 
            call. = FALSE
        )
    }
    
    # check selected indicators
    valid_indicators <- c("N", "LNN", "S", "D", "SN", "SNA", "H", 
                          "AMBI", "ITI", "L", "PIE", "N2")
    # total abundance, natural log of total abundance,
    # species sensitivity, Margalef's D, Rygg's SN and its adjustment SNA, 
    # Shannon H, AMBI, ITI, Simpson's L, Hurlbert's PIE, Hill's N2
    if (length(settings$indicators) == 0L) {
        stop(
            sprintf(
                fmt = paste0(
                    "No indices have been specified in JSON-file %s.\n",
                    "Please choose a subset of %s (see package vignette)"),
                sQuote(basename(filename)),
                toString(valid_indicators)
            ), 
            call. = FALSE
        )
    }
    if (length(settings$indicators) > 3L) {
        stop(
            sprintf(
                fmt = "Maximum 3 indices allowed (%i found in %s).",
                length(settings$indicators),
                sQuote(basename(filename))
            ), 
            call. = FALSE
        )
    }
    is_invalid <- !(tolower(settings$indicators) %in% tolower(valid_indicators))
    if (any(is_invalid)) {
        stop(
            sprintf(
                fmt = paste0(
                    "Invalid indices found in JSON-file %s: %s\n",
                    "Please choose a subset of %s (see package vignette)"),
                sQuote(basename(filename)),
                toString(sQuote(settings$indicators[is_invalid])), 
                toString(valid_indicators)
            ), 
            call. = FALSE
        )
    }
    settings$indicators <- tolower(settings$indicators)

    # check indicator weights
    if (!is.null(settings$weights)) {
        w <- settings$weights
        if (is.character(w)) { 
            # handle rational numbers like 1/3 etc.
            tmp <- w
            w <- sapply(X = parse(text = w), FUN = eval)
            attr(w, "character") <- tmp
        }
        if (length(w) != length(settings$indicators)) {
            stop(
                sprintf(
                    fmt = paste("The number of weights (%i) should be equal",
                        "to the number of indices (%i)"),
                    length(w),
                    length(settings$indicators)                
                ), 
                call. = FALSE
            )
        }
        is_negative <- w < -.Machine$double.eps
        if (any(is_negative)) {
            stop(
                sprintf(
                    fmt = "Negative weights found in JSON-file: %s",
                    sQuote(basename(filename))
                ), 
                call. = FALSE
            )
        }
        if (abs(1 - sum(w)) > 1.0e-3) {
            message("Weights do not sum to one and will be normalized first.")
            w <- abs(w) # small neg.numbers -> 0
            w <- w / sum(w)
            attr(w, "normalized") <- TRUE
        } else {
            attr(w, "normalized") <- FALSE
        }
        settings$weights <- w
    }
    
    # check confidence level
    is_invalid <- 
        (settings$confidencelevel < 0.50) ||
        (settings$confidencelevel > 0.99)
    if (is_invalid) {
        stop(
            sprintf(
                fmt = "Confidence level needs to be in [0.5, 0.99]. Check JSON-file: %s",
                sQuote(basename(filename))
            ), 
            call. = FALSE
        )
    }

    # process pressure info
    if (is.null(settings$pressure)) {
        if (is.null(settings$weights)) {
            stop(
                sprintf(
                    fmt = "weights should be given in %s when pressure is missing", 
                    sQuote(basename(filename))
                ), 
                call. = FALSE
            )
        }
    } else {
        required_keys <- c("name", "unit")
        names(settings$pressure) <- tolower(names(settings$pressure))
        found <- tolower(required_keys) %in% names(settings$pressure)
        if (any(!found)) {
            stop(
                sprintf(
                    fmt = "key pressure.%s is missing in %s\n(see package vignette)", 
                    toString(sQuote(required_keys[!found])), 
                    sQuote(basename(filename))
                ), 
                call. = FALSE
            )
        }
    }
    

    # check months
    if (!is.integer(settings$months) || length(settings$months) != 2L) {
    	stop(
            "key 'months' should be an integer vector of length 2", 
             call. = TRUE
        )
    }
	if (!all(settings$months %in% 1:12)) {
		stop("elements of key 'months' should be in [1, 12]", call. = TRUE)
	}
    if ((settings$months[2] - settings$months[1]) < 0) {
		stop(
            "First month to analyse should be smaller than or equal to last month", 
            call. = TRUE
        )
	}

    # check data pooling
    names(settings$pooling) <- tolower(names(settings$pooling))
    if (!is.logical(settings$pooling$enabled)) {
        stop(
            "key 'pooling:enabled' should be either 'true' or 'false'", 
            call. = TRUE
        )
    }
    if (settings$pooling$enabled) {
        settings$pooling$targetarea <- range(settings$pooling$targetarea)
        if (!is.numeric(settings$pooling$targetarea) || 
            (length(settings$pooling$targetarea) != 2L)) {
        	stop(
                "key 'pooling:targetArea' should be a numeric vector of length 2", 
                call. = TRUE
            )
        }
        settings$pooling$randomseed <- as.integer(settings$pooling$randomseed)
        if (!is.integer(settings$pooling$randomseed)) {
            stop(
                "key 'pooling:randomSeed' needs to be an integer vector of length 1", 
                call. = TRUE
            )
        }
    }

    # genus to species conversion
    if (!is.logical(settings$genustospeciesconversion)) {
        stop(
            "key 'genusToSpeciesConversion' should be either 'true' or 'false'", 
            call. = TRUE
        )
    }
    
    # check key 'legendText'
    if (is.null(settings$legendtext)) {
        settings$legendtext <- "normalized"
    } else {
        if ( (length(settings$legendtext) != 1L) ||
            !(tolower(settings$legendtext) %in% c("eqr", "normalized"))) {
            stop(
                sprintf(
                    "key %s should be either %s or %s",
                    sQuote("legendText"),
                    sQuote("EQR"),
                    sQuote("normalized")
                ),
                call. = FALSE
            )
        }
    }

    # check key 'model'
    if (is.null(settings$model)) {
        settings$model <- "linear"
    } else {
        if ( (length(settings$model) != 1L) ||
            !(settings$model %in% c("linear", "exponential"))) {
            stop(
                sprintf(
                    "key %s should be either %s or %s",
                    sQuote("model"),
                    sQuote("linear"),
                    sQuote("exponential")
                ),
                call. = FALSE
            )
        }
    }

    
    # issue a warning if keys in JSON-file have not been used
    not_used <- setdiff(tolower(names(settings)), tolower(available_keys))
    if (length(not_used) > 0L) {
        warning(
            sprintf(
                fmt = "key(s) %s not used in JSON-file: %s",
                toString(not_used), 
                sQuote(basename(filename))
            ),
            call. = FALSE
        )
    }
    
	# return results
	settings
}



#'  Read and Validate BENMMI Input Files
#'
#' 	This function reads and checks benthos files. The format is a superset
#'  of the BEQI2-format as specified in Van Loon (2013). In addition to the
#'  BEQI2-format, the benthos-format also includes columns latitude (LAT), 
#'  longitude (LONG), and sieve mesh size (MESH).
#'
#' @param filename name of benthos file (\code{character})
#' 
#' @import benthos
#' @importFrom purrr walk
#'
#' @references Willem van Loon, 2013. BEQI2 INPUT FORMAT
#'  
#' @seealso \code{\link[benthos]{read_beqi2}}
read_mmi <-
function(filename) {

    # read BEQI2-format
    d <- read_beqi2(filename)

    # additional validation
    required_vars <- c("LAT", "LONG", "MESH")
    missing_vars <- setdiff(required_vars, toupper(names(d)))
	if (length(missing_vars) > 0L) {
		stop(
            sprintf(
                fmt = "The following columns are missing: %s", 
                toString(missing_vars)
            ),
            call. = FALSE
        )
	}
    
    # check availability of pressure column
    if ("PRESSURE" %in% toupper(names(d))) {
        required_vars <- c(required_vars, "PRESSURE")
        
        #  check values of pressure column
        if (!is.numeric(d$PRESSURE)) {
            d$PRESSURE <- try(
                as.Date(d$PRESSURE), 
                silent = TRUE
            )
            if (inherits(d$PRESSURE, "try-error")) {
                stop(
                    "Either dates (YYYY-MM-DD) or numeric values expected\n",
                    "in the PRESSURE column of the benthos file", 
                    call. = FALSE
                )
            }
            # convert date to numeric value
            year <- d$PRESSURE  %>% 
                format("%Y") %>% 
                as.numeric
            day_of_year <- d$PRESSURE %>% 
                format("%j") %>% 
                as.numeric
            days_in_year <- paste0(year, "-12-31") %>% 
                as.Date %>% 
                format("%j") %>% 
                as.numeric
            d$PRESSURE <- year + (day_of_year - 0.5) / days_in_year
        }
    }
    
    # check completeness of required variables
    required_vars %>%
        walk(function(x) {
            is_na <-is.na(d[[x]])
            if (any(is_na)) {
        		stop(
                    sprintf(
                        fmt = "Column %s in %s contains missing values at lines: %s",
                        sQuote(x),
                        sQuote(filename),
                        toString(which(is_na))
                    ),
                    call. = FALSE
                )
            }
        })

    # check uniqueness of coordinates
    n1 <- d %>% 
        select_(~OBJECTID, ~SAMPLEID, ~DATE) %>% 
        distinct_ %>% 
        nrow
    n2 <- d %>% 
        select_(~OBJECTID, ~SAMPLEID, ~DATE, ~LAT, ~LONG) %>% 
        distinct_ %>% 
        nrow
    if (n1 != n2) {
		stop(
            "Coordinates have to be unique for combinations of ",
            "OBJECTID, SAMPLEID, DATE",
            call. = FALSE
        )
	} 

    # check uniqueness of sieve mesh
    mesh_size <- unique(d$MESH)
    is_unique <- length(mesh_size) == 1L
    if (!is_unique) {
		stop(
            sprintf(
                fmt = "Sieve mesh size is not unique.\nThe following mesh sizes have been found in column MESH: %s", 
                toString(mesh_size)
            ),
            call. = FALSE
        )
    }

    # return result
    d
}



#' Perform BENMMI Analysis
#'
#' This function performs a complete BENMMI analysis following the
#' settings provided in \code{filename}.
#'
#' @param filename name of the JSON file defining all analysis steps.
#' @param tmpdir directory to store temporary files (for debugging only)
#' @param browse load resulting report in a browser? \code{TRUE} or \code{FALSE}
#'
#' @import xtable
#' @import dplyr
#' @import markdown
#' @importFrom tidyr spread gather extract_numeric unnest
#' @import knitr
#' @import tcltk
#' @import ggplot2
#' @importFrom utils browseURL flush.console
#' @importFrom readr read_csv cols_only col_character
#' 
#' @examples
#' 
#'\donttest{
#' # This example illustrates a typical use case of the BENMMI-package.
#' # Note: execution may take several minutes.
#' # See the package vignette for more advanced examples and details. 
#' 
#' if (interactive()) {
#' 
#' # Create a work directory (in this example, a temporary
#' # directory, but in real use cases a persistent directory 
#' # will obviously be more useful).
#' my_dir <- tempfile("benmmi-example")
#' dir.create(my_dir)
#' 
#' # Populate this directory with simple use cases 
#' # (see the package-vignette for details).
#' # Most users will probably use one of these use cases as a 
#' # template for their own study.
#' BENMMIdir(my_dir)
#' 
#' # Run BENMMI given the settings in "settings-S-D-lin.json". This file
#' # relates to one of the predefined use cases. 
#' my_settings_file <- file.path(my_dir, "settings-S-D-lin.json")
#' benmmi(my_settings_file, browse = FALSE)
#' 
#' # The output (HTML-report and data-files) is stored in 'my_dir'
#' # and described in the package-vignette and resulting HTML-report itself.
#' # It is also possible to directly view the generated
#' # HTML-report by setting the browse-argument of the benmmi-function to TRUE.
#' }
#' }
#' @rdname benmmi-main
#' @export
benmmi <-
function(filename = NULL, tmpdir = tempfile(pattern = "BENMMI"), browse = TRUE) {

    # prevent potential problems with dates in other locales
    old_locale <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_TIME", old_locale))
    Sys.setlocale("LC_TIME", "C")
    
    # interactive selection of filename
    if (is.null(filename)) {
        if (capabilities("tcltk")) {
            filename <- tk_choose.files(
                default = "", 
                caption = "Select file with BENMMI settings",
                multi = FALSE, 
                filters = matrix(data = c("BENMMI settings", ".json"), nrow = 1)
            )
        } else {
            stop(
                "The 'tcltk'-package is not supported on this machine.\n",
                "Please provide a valid filename as function argument\n",
                call. = FALSE
            )
        }
    }

    # stop if the user presses Cancel or Esc
    if(length(filename) == 0L) {
        message("The BENMMI run has been cancelled by the user.")
        return(invisible(NULL))
    }
    
    # check if filename exists
    if (!file.exists(filename)) {
        stop(
            sprintf("JSON-file %s does not exist", sQuote(filename)), 
            call. = FALSE
        )
    }

    # initialization message
    message("The BENMMI tool is running...")
    flush.console()

    # read settings
	settings <- read_settings(filename)

    # set working directory
    owd <- getwd()
    on.exit(setwd(owd), add = TRUE)
    setwd(dirname(filename))

    # normalize paths (full paths to make package more robust)
    for (f in names(settings$files)) {
        settings$files[[f]] <- suppressWarnings(normalizePath(settings$files[[f]]))
    }

    # add output files
    output_dir <- file.path(
        getwd(), 
        paste0("OUTPUT-", format(Sys.time(), format = "%Y%m%dT%H%M%S"))
    )
    dir.create(output_dir)
    prefix <- sub(
        pattern = "\\.[[:alnum:]]+$", 
        replacement = "", 
        x = basename(settings$files$benthos)
    )
    settings$files$log <- file.path(output_dir, 
        paste0("LOG-", prefix, ".log"))
    settings$files$out_sample <- file.path(output_dir, 
        paste0("SAMPLE-", prefix, ".csv"))
    settings$files$out_habitat <- file.path(output_dir, 
        paste0("HABITAT-", prefix, ".csv"))
    settings$files$out_objectid <- file.path(output_dir, 
        paste0("OBJECTID-", prefix, ".csv"))
    settings$files$out_group <- file.path(output_dir, 
        paste0("GROUP-", prefix, ".csv"))
    settings$files$out_iti <- file.path(output_dir, 
        paste0("ITI-", prefix, ".csv"))
    settings$files$out_tidy <- file.path(output_dir, 
        paste0("TIDY-", prefix, ".csv"))
    if (settings$pooling$enabled) {
        settings$files$pooling <- file.path(output_dir, 
            paste0("POOLING-", prefix, ".csv"))
    }
    settings$files$report <- file.path(output_dir, 
        paste0("REPORT-", prefix, ".html"))
    # settings$files$snap_shot <- file.path(output_dir,
    #     paste0(sub(pattern = "OUTPUT", replacement = "BENMMI-snapshot", 
    #                x = basename(output_dir)), ".zip"))

    # start logging
    to_log <- function(level = c("INFO", "WARNING", "ERROR"), message) {
        level <- match.arg(level)
        cat(
            format(Sys.time()), " [", level, "] ", message, "\n", 
            sep = "",
            file = settings$files$log, 
            append = TRUE
        )
        if (level != "INFO") {
            switch(level,
               "ERROR"   = stop(message, call. = FALSE),
               "WARNING" = warning(message, call. = FALSE)
            )
        }
    }
    to_log("INFO", "Starting a new BENMMI session")
    on.exit(to_log("INFO", "The BENMMI session has been terminated"), add = TRUE)

    # initialize random number generator
    if (settings$pooling$enabled) {
        to_log("INFO", "Initializing the pseudo random number generator...")
        if (is.null(settings$pooling$randomseed)) {
            to_log("INFO", "No seed has been specified.")
            to_log("INFO", "The default initialization process will be followed.")
        } else {
            set.seed(seed = settings$pooling$randomseed)
        }
        to_log("INFO", "the pseudo random number generator has been initialized.")
    }

    # check existence of taxa-file
    to_log("INFO", 
          sprintf(
                "Checking the existence of taxa-file %s...",
                sQuote(basename(settings$files$taxa))
          )
    )
    if (!file.exists(settings$files$taxa)) {
        to_log("ERROR", "The taxa-file has not been found")
        return(invisible(NULL))
    }
    to_log("INFO", "the taxa-file has been found")

    # read taxa-file
    to_log("INFO", "Reading the taxa-file...")
    d_taxa <- tryCatch(
        read_taxa(filename = settings$files$taxa),
        error = function(e) {
            to_log("ERROR", sprintf("while reading taxa-file. %s", e$message))
        }
    )
    to_log("INFO", "the taxa-file has been read")

    # check existence of 'groupsToExclude'-file
    to_log("INFO", 
          sprintf(
                "Checking the existence of 'groups to exclude'-file %s...",
                sQuote(basename(settings$files$groupstoexclude))
          )
    )
    if (!file.exists(settings$files$groupstoexclude)) {
        to_log("ERROR", "The 'groupsToExclude'-file has not been found")
        return(invisible(NULL))
    }
    to_log("INFO", "the 'groupsToExclude'-file has been found")

    # read 'groupsToExclude'-file
    to_log("INFO", "Reading the 'groupsToExclude'-file...")
    d_groups <- tryCatch(
        settings$files$groupstoexclude %>%
            read_csv(
                col_types = cols_only(
                    GROUP = col_character(),
                    DESCRIPTION = col_character()
                )
            ),
        error = function(e) {
            to_log("ERROR", sprintf("while reading taxa-file. %s", e$message))
        }
    )
    to_log("INFO", "the 'groupsToExclude'-file has been read")

    # check existence of the benthos file
    to_log("INFO", 
          sprintf(
                "Checking the existence of benthos file %s...",
                sQuote(basename(settings$files$benthos))
          )
    )
    if (!file.exists(settings$files$benthos)) {
        to_log("ERROR", "The benthos file has not been found")
        return(invisible(NULL))
    }
    to_log("INFO", "the benthos file has been found")

    # read benthos-file
    to_log("INFO", "Reading the benthos-file...")
    d_mmi <- tryCatch(
        read_mmi(filename = settings$files$benthos),
        error = function(e) {
            to_log("ERROR", sprintf("while reading benthos-file. %s", e$message))
        }
    )
    if (!("PRESSURE" %in% names(d_mmi)) || all(is.na(d_mmi$PRESSURE))) {
        d_mmi$PRESSURE <- NULL
        if (is.null(settings$weights)) {
            stop(
                sprintf(
                    fmt = "weights should be given in %s when PRESSURE column is missing in %s", 
                    sQuote(basename(filename)),
                    sQuote(basename(settings$files$benthos))
                ), 
                call. = FALSE
            )
        }
    }
    if (!("PRESSURE" %in% names(d_mmi))) {
        if (!is.null(settings$pressure)) {
            to_log("WARNING",
                sprintf(
                    fmt = paste(
                        "pressure has been specified in %s but is missing in %s\n",
                        "BENMMI continues without pressure optimization"),
                    sQuote(basename(filename)),
                    sQuote(basename(settings$files$benthos))
                )
            )
            settings$pressure <- NULL
        }
    }
    to_log("INFO", "the benthos file has been read")
    
    # check if records are within the period of interest
    in_poi <- d_mmi$DATE %>% 
        format(format = "%m") %>% 
        as.integer %>% 
        between(settings$months[1], settings$months[2])
    if (!any(in_poi)) {
        to_log("ERROR", 
            sprintf("No months in file %s are in the specified interval [%s].",
                sQuote(basename(settings$files$benthos)),
                paste(settings$months, collapse = ", ")
            )
        )
    }

    # make sure that names in benthos file correspond to those in the taxa-list
    d_mmi <- d_mmi %>% 
        mutate_(
            TAXON_OLD = ~TAXON, 
            TAXON = ~as_accepted(taxon = TAXON, taxa = d_taxa)
        ) 
    
    # add unique sampling unit identifier
    d_mmi <- d_mmi %>%
        select_(~OBJECTID, ~SAMPLEID, ~DATE) %>%
        distinct_ %>%
        mutate_(ID = ~row_number()) %>%
        inner_join(d_mmi, by = c("OBJECTID", "SAMPLEID", "DATE"))

    # conditionally read user-defined AMBI
    d_ambi <- NULL
    if (("ambi" %in% settings$indicators) &&
        !is.null(settings$files$ambi) && 
        (settings$files$ambi != "")) {
        to_log("INFO", 
              sprintf(
                    "Checking the existence of AMBI-file %s...",
                    settings$files$ambi %>% basename %>% sQuote
              )
        )
        if (!file.exists(settings$files$ambi)) {
            to_log("ERROR", "the AMBI-file has not been found")
            return(invisible(NULL))
        }
        to_log("INFO", "the AMBI-file has been found")
        to_log("INFO", "Reading the AMBI-file...")
        d_ambi <- tryCatch(
            read_ambi(filename = settings$files$ambi),
            error = function(e) {
                to_log("ERROR", sprintf("while reading AMBI-file. %s", e$message))
            }
        )
        to_log("INFO", "the AMBI-file has been read")
        
        # make sure that names in sensitivity file correspond to those in the taxa-list
        d_ambi <- d_ambi %>% 
            mutate_(TAXON = ~as_accepted(taxon = TAXON, taxa = d_taxa)) %>%
            distinct_
        
        # check if taxa are still unique after conversion to WoRMS
        if (anyDuplicated(d_ambi$TAXON)) {
            to_log("WARNING",
                paste0(
                    "the AMBI-file causes inconsistencies\n",
                    "(TAXON-AMBI class combinations are not unique)\n",
                    "Only the first combination will be used"
                )
            )
        }
        d_ambi <- d_ambi %>% distinct_(.dots = "TAXON", .keep_all = TRUE)
    }

    # read user-defined ITI
    if ("iti" %in% settings$indicators) {
        if (is.null(settings$files$iti) || (settings$files$iti == "")) {
            to_log("ERROR", "the ITI-file has not been specified")
        }
        to_log("INFO", 
              sprintf(
                    "Checking the existence of ITI-file %s...",
                    settings$files$iti %>% basename %>% sQuote
              )
        )
        if (!file.exists(settings$files$iti)) {
            to_log("ERROR", "the ITI-file has not been found")
        }
        to_log("INFO", "the ITI-file has been found")
        to_log("INFO", "Reading the ITI-file...")
        d_iti <- tryCatch(
            read_iti(filename = settings$files$iti),
            error = function(e) {
                to_log("ERROR", sprintf("while reading ITI-file. %s", e$message))
            }
        )
        to_log("INFO", "the ITI-file has been read")
        
        # make sure that names in sensitivity file correspond to those in the taxa-list
        d_iti <- d_iti %>% 
            mutate_(TAXON = ~as_accepted(taxon = TAXON, taxa = d_taxa)) %>%
            distinct_
        
        # check if taxa are still unique after conversion to WoRMS
        if (anyDuplicated(d_iti$TAXON)) {
            to_log("WARNING",
                paste0(
                    "the ITI-file causes inconsistencies\n",
                    "(TAXON-ITI class combinations are not unique)\n",
                    "Only the first combination will be used"
                )
            )
        }
        d_iti <- d_iti %>% distinct_(.dots = "TAXON", .keep_all = TRUE)
    }
    
    # read habitat reference file
    to_log("INFO", 
          sprintf(
                "Checking the existence of habitat reference file %s...",
                settings$files$habitats %>% basename %>% sQuote
          )
    )
    if (!file.exists(settings$files$habitats)) {
        to_log("ERROR", "the habitat reference file has not been found")
        return(invisible(NULL))
    }
    to_log("INFO", "the habitat reference file has been found")
    to_log("INFO", "Reading the habitat reference file...")
    d_ref <- tryCatch(
        read_ref(filename = settings$files$habitats, settings$indicators),
        error = function(e) {
            to_log("ERROR", 
                sprintf(
                    paste0(
                        "while reading habitat reference file:\n%s\n",
                        "see vignette for its format and how to estimate reference values."
                    ),
                    e$message
                )
            )
        }
    )
    to_log("INFO", "the habitat reference file has been read")

    # check if reference data are available for all records in d_mmi
    to_log(
        "INFO", 
        "Checking if reference data are available for all records in the benthos-file..."
    )
    d <- d_mmi %>% 
        select_(~OBJECTID, ~HABITAT) %>% 
        distinct_ %>%
        left_join(d_ref, by = c("OBJECTID", "HABITAT"))
    col_names <- names(d) %>% 
        setdiff(c("OBJECTID", "HABITAT"))
    d <- d %>% 
        filter_(.dots = paste(sprintf("is.na(%s)", col_names), collapse = "|"))
    if (nrow(d) > 0L) {
        M <- NULL
        for (i in 1:nrow(d)) {
            m <- d[i, ] %>% 
                unlist(use.names = TRUE)
            m2 <- m[is.na(m)]
            m1 <- m[c("OBJECTID", "HABITAT")]
            m1 <- paste(names(m1), sQuote(m1), sep = "=", collapse = " AND ") 
            m <- paste0(paste(names(m2), collapse = ", "), " for ", m1)
            M <- c(M, m)
        }
        to_log("ERROR", sprintf(
                "The following columns are empty in the habitat reference file:\n%s",
                paste(M, collapse = ";\n")
            )
        )
        return(invisible(NULL))
    }
    to_log("INFO", "reference data are available for all records in the benthos-file.")
    
	# create temporary directory
	if (!file.exists(tmpdir)) {
        to_log("INFO", "Creating a temporary directory...")
		dir.create(tmpdir)
	}
    to_log("INFO", "a temporary directory has been created.")

	# copy template of the report to temporary directory
    to_log("INFO", "Populating the temporary directory...")
	templates <- list.files(
		path = system.file("Rmd", package = "BENMMI"),
		pattern = "\\.Rmd$", full.names = TRUE)
	file.copy(from = templates, to = tmpdir)
    to_log("INFO", "the temporary directory has been populated.")

    # create Markdown document 
    # (code below works better than knit2html)
    to_log("INFO", "Starting to create a report...")
    setwd(tmpdir)
	suppressMessages(
        mdfile <- try(knit(input = "benmmi.Rmd", quiet = TRUE), silent = TRUE)
    )
    if (inherits(mdfile, "try-error")) {
        to_log(
            level = "ERROR", 
            message = toString(attr(mdfile, "condition")$message)
        )
        return(invisible(NULL))
    }
    to_log("INFO", "a report has been created.")
    to_log("INFO", "Converting the report to HTML...")
	output <- markdownToHTML(
		file  = mdfile, 
		output = NULL,
        options = getOption("markdown.HTML.options"),
        extensions = getOption("markdown.extensions"),
    	title = "BENMMI Report",
        stylesheet = system.file("css", "benmmi.css", package = "BENMMI")
	)
    writeLines(text = output, con = settings$files$report)
    to_log("INFO", "the report has been converted to HTML.")
    
    # create data snap shot of outputs
    # zip(
    #     zipfile = settings$files$snap_shot, 
    #     files = c(settings$files$log, settings$files$out_sample, 
    #               settings$files$out_habitat, settings$files$out_objectid, 
    #               settings$files$report, settings$files$pooling), 
    #     flags = "-qj9X"
    # )

	# view result
	if (browse) {
		browseURL(settings$files$report)
	}
    
    # finalization
    message("The BENMMI run has been completed successfully.")
}



#' Perform BENMMI Analysis
#'
#' @rdname benmmi-main
#'  
#' @export
BENMMI <- 
function(filename = NULL, tmpdir = tempdir(), browse = TRUE) {
    benmmi(filename = filename, tmpdir = tmpdir, browse = browse)
}



#'  Create BENMMI Directory Structure
#'
#'  Creates a BENMMI-directory structure and populates it with some
#'  relevant BENMMI-files. Users may wish to modify this directory structure 
#'  and add their own data.
#'
#' @param path name of an exisiting directory. This directory should
#'      be empty to prevent loss of data. If missing, a dialogue will
#'      appear.
#'  
#' @export
BENMMIdir <- function(path = NULL) {

    # interactive selection
    if (is.null(path)) {
        if (capabilities("tcltk")) {
            path <- tk_choose.dir(
                caption = "Select directory to store BENMMI files"
            )
        } else {
            stop(
                "The 'tcltk'-package is not supported on this machine.\n",
                "Please provide a valid path as function argument\n",
                call. = FALSE
            )
        }
    }
    
    # check path
    if (is.na(path)) {
        message("The BENMMI run has been cancelled by the user.")
        return(invisible(NULL))
    }
    if (!file.exists(path)) {
        stop("directory does not exist", call. = FALSE)
    }
    
    # check if directory is empty (to prevent overwriting existing files)
    if (length(list.files(path)) != 0L) {
        stop(sprintf("directory %s is not empty!", sQuote(path)), call. = FALSE)
    }

    # populate directories
    tmp <- file.copy(
        from = c(
            system.file("extdata/INPUT-FILES", package = "BENMMI"),
            system.file("extdata/REF-FILES", package = "BENMMI"),
            system.file("extdata/settings-S-D-lin.json", package = "BENMMI"),
            system.file("extdata/settings-D-exp.json", package = "BENMMI")
        ),
        to = path,
        recursive = TRUE
    )

    # show message
    message(
        sprintf(
            paste0(
                "Directory %s\nhas been populated with BENMMI files.\n",
                "To run the BENMMI tool, type: BENMMI() or benmmi()\n",
                'For the tutorial, type: vignette("BENMMI_User_Manual")\n',
                "For more technical information, type: ?benmmi"
            ), 
            sQuote(path)
        )
    )
}

