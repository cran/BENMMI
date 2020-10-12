## ----ini, echo=FALSE, results='hide', message=FALSE---------------------------
library(BENMMI)
library(benthos)
library(knitr)
library(xtable)
library(ggplot2)
library(DEoptim)
library(readr)
library(dplyr)

## ----echo=FALSE---------------------------------------------------------------
opts_chunk$set(
    echo = FALSE,
    comment = NA,
    quiet = TRUE,
    progress = FALSE,
    tidy = FALSE,
    cache = FALSE,
    message = FALSE,
    error = TRUE,
    warning = TRUE
)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  library(BENMMI)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  BENMMIdir()

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  BENMMIdir(path = "c:/myprojects/BENMMI/BENMMI_FILES")

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  BENMMI()

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  BENMMI(filename = "c:/myprojects/BENMMI/BENMMI_FILES/settings-S-D-lin.json")

## ----echo=FALSE, results='asis'-----------------------------------------------
cat(
    paste(
        readLines(system.file("extdata", "settings-S-D-lin.json", package = "BENMMI")), 
        collapse = "\n"
    )
)

## ----echo=FALSE, results='asis'-----------------------------------------------
d <- scan(
    file = "./tables/tab-benmmi-input.csv", 
    what = character(), 
    sep = ",",
    quiet = TRUE
)
h <- d[1:3]
d <- as.data.frame(matrix(data = d[-(1:3)], ncol = 3, byrow = TRUE))
colnames(d) <- h
print(
    xtable(x = d, align = "llp{50mm}p{50mm}"), 
    include.rownames = FALSE, 
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----echo=FALSE, results='asis'-----------------------------------------------
d <- read_ambi(
    filename = system.file(
        "extdata", "REF-FILES", "AMBI-NL+.csv", 
        package = "BENMMI"
    )
)
print(
    xtable(x = d[sample.int(n = nrow(d), size = 25), ], align = "llr"), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----echo=FALSE, results='asis'-----------------------------------------------
d <- read_iti(
    filename = system.file(
        "extdata", "REF-FILES", "ITI+carnivores-2015-10-23.csv", 
        package = "BENMMI"
    )
)
print(
    xtable(x = d[sample.int(n = nrow(d), size = 25), ], align = "llrr"), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----echo=FALSE, results='asis'-----------------------------------------------
filename <- system.file(
    "extdata", "REF-FILES", "AREAS-HABITATS-SNS-2016-11-27.csv", 
    package = "BENMMI"
)
d <- read_ref(file = filename, indicators = c("D", "S", "AMBI", "ITI"))
print(
    xtable(x = d), 
    include.rownames = FALSE,
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}"),
    sanitize.text.function = function(x) {
        gsub(pattern="_", replacement = "\\\\_", x=x)
    },
    rotate.colnames = TRUE
)

## ----echo=FALSE, results='asis'-----------------------------------------------
filename <- system.file(
    "extdata", "REF-FILES", "TAXONOMIC-GROUPS-EXCLUDED.csv", 
    package = "BENMMI"
)
d <- read_csv(
    file = filename, 
    col_types = cols(
      GROUP = col_character(),
      DESCRIPTION = col_character()
    ))
print(
    xtable(x = d, align = "lll"), 
    include.rownames = FALSE, 
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----echo=FALSE, results='asis'-----------------------------------------------
filename <- system.file(
    "extdata", "REF-FILES", "TAXA-BE-DE-NL-UK-2017-01-06.csv", 
    package = "BENMMI"
)
d <- read_csv(
    file = filename, 
    col_types = cols(
      group = col_character(),
      provided = col_character(),
      accepted = col_character(),
      level = col_character(),
      quality_code = col_integer()
    )
)
print(
    xtable(x = d %>% slice(1:10), align = "lllllr"), 
    include.rownames = FALSE, 
    size = "footnotesize",
    add.to.row = list(list(-1), "\\rowcolor{blue!15}")
)

## ----echo=FALSE, results='asis'-----------------------------------------------
# https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
d <- read.csv(file = "./tables/tab-GROUP.csv")
print(xtable(d[1:20, 1:11]), floating.environment = 'sidewaystable', size = 'footnotesize', include.rownames = FALSE)

## ----echo=FALSE, results='asis'-----------------------------------------------
# https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
d <- read.csv(file = "./tables/tab-ITI.csv")
names(d)[ncol(d)] <- "NA"
print(xtable(d[1:20, ]), floating.environment = 'sidewaystable', size = 'footnotesize', include.rownames = FALSE)

## ----echo=FALSE, message=FALSE------------------------------------------------
set.seed(314)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
#  #Calculated with BENMMI 2016-08-18
#  d <- read_csv("saltkallefjord-indices.csv") %>%
#      mutate(
#          date = SAMPLEID %>%
#              substr(start = 7, stop = 16) %>%
#              as.Date
#      ) %>%
#      select(date, D)
#  d %>%
#      arrange(date) %>%
#      write_csv("saltkallefjord-D.csv")

## ----echo=FALSE---------------------------------------------------------------

plot_fit <- function(d, theta = NULL) {
    g <- ggplot() + 
        geom_point(
            data = d,
            mapping = aes(x = date, y = D), 
            colour = "blue"
        ) +
        scale_x_date(
            name = "", 
            limits = as.Date(c("1965-01-01", "1980-01-01"))
        ) +
        scale_y_continuous(name = "Margalef diversity D")
    if (is.null(theta)) {
        return(g)
    }
    d_fit <- data.frame(
        date = seq(
            from = min(d$date), 
            to = max(d$date), 
            by = 0.1
        )
    )
    d_fit$ndate <- as.numeric(d_fit$date)
    d_fit$D <- f_gl(x = d_fit$ndate, theta)
    g +
        geom_path(
            data = d_fit, 
            mapping = aes(x = date, y = D)
        )
}

## ----echo=TRUE----------------------------------------------------------------
f_gl <- function(x, theta) {
    A <- theta[1]
    K <- theta[2]
    B <- theta[3]
    M <- theta[4]
    A + (K - A) / (1 + exp(-B*(x-M)))
}

## ----echo=TRUE----------------------------------------------------------------
# read data
d_obs <- read.csv("./data/saltkallefjord-D.csv", as.is = TRUE)

# coercion from character to Date-object
d_obs$date <- as.Date(d_obs$date)

# print head of these data
head(d_obs)

## ----fig.width=4, fig.height=3, out.width="0.7\\textwidth", echo=FALSE--------
plot_fit(d_obs)

## ----echo=TRUE----------------------------------------------------------------
d_obs$ndate <- as.numeric(d_obs$date)

## ----echo=TRUE----------------------------------------------------------------
f_obj <- function(theta, data) {
    # constraint: lower asymptote should be nonnegative
    if (theta[1] < 0) {
        return(Inf)
    }
    
    # predict Margalef diversity by means of the generalised logistic function
    D_hat <- f_gl(data$ndate, theta)
    
    # difference between Margalef diversity based on 
    # observations and the generalised logistic function
    error <- d_obs$D - D_hat
    
    # our objective to minimize: the sum of squared errors
    sum(error * error)
}

## ----echo=TRUE----------------------------------------------------------------
# define lower and upper bounds for each parameter
lower_bounds <- c(0, 5, 0.0001, -1000)
upper_bounds <- c(5, 9, 1.0000,  1000)

# DE-optimization
opt <- DEoptim(
    fn = f_obj, 
    lower = lower_bounds, 
    upper = upper_bounds, 
    data = d_obs, 
    control = DEoptim.control(NP = 100, itermax = 100, trace=FALSE)
)

## ----echo=FALSE---------------------------------------------------------------
# check results
theta <- opt$optim$bestmem
in_bounds <- all(abs(theta - lower_bounds) > 1.0e-3) & 
             all(abs(theta - upper_bounds) > 1.0e-3)
    
if (!in_bounds) {
    cat("parameters pressed against bounds: ", toString(p), "\n")
    print(lower_bounds)
    print(upper_bounds)
}


## ----fig.width=4, fig.height=3, out.width="0.7\\textwidth", echo=FALSE--------
plot_fit(d_obs, theta)

## ----echo=FALSE, results='asis'-----------------------------------------------
# https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf
d <- read.csv(file = "./tables/tab-checks.csv")
print(
    xtable(d, align = "llp{55mm}lp{30mm}"), 
    hline.after = c(-1), 
    add.to.row = list(
        pos = list(0),
        command = paste0(
            "\\hline\n\\endhead\n",
            "\\hline\n",
            "\\multicolumn{", ncol(d) + 1, "}{l}",
            "{\\footnotesize Continued on next page}\n",
            "\\endfoot\n",
            "\\endlastfoot\n")),
    floating = FALSE,
    tabular.environment = "longtable"
)


