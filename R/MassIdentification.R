#' @title MassIdentification
#' @description Identify feature in HMDB according to mz.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param mass.tolerance mz tolerance.
#' @param polarity polarity
#' @return A MetFlowData added HMDB result in tags.
#' @export
#' @details \href{http://www.hmdb.ca/}{HMDB} is a metabolite database. The
#' database data can be loaded. In positive mode, adductions of +H, +NH4,
#'  +Na, +H-2H2O, +H-H2O, +K have been considered, and in negative mode, -H,
#'  +Cl, +CH3COO, -H-H2O, +F have been considered.

MassIdentification <- function(MetFlowData = MetFlowData,
                               mass.tolerance = 30,
                               polarity = "positive",
                               show = 5) {
  tags <- MetFlowData[["tags"]]
  mz1 <- as.numeric(tags[, "mz"])
  data("hmdbdatabase")
  mz2 <- hmdbdatabase[, "Mass"]
  HMDB.ID <- hmdbdatabase[, 1]
  HMDB.name <- hmdbdatabase[, 2]
  HMDB.formula <- hmdbdatabase[, 4]

  if (polarity == "positive") {
    mz2.h <- mz2 + 1.0078
    mz2.nh4 <- mz2 + 18.0344
    mz2.na <- mz2 + 22.9898
    mz2.h_2h2o <- mz2 + 1.0078 - 2 * 18.0106
    mz2.h_h2o <- mz2 + 1.0078 - 18.0106
    mz2.k <- mz2 + 38.9637
    mz2.all.pos <-
      data.frame(mz2.h, mz2.nh4, mz2.na, mz2.h_2h2o, mz2.h_h2o, mz2.k)
  }

  if (polarity == "negative") {
    mz2_h <- mz2 - 1.0078
    mz2.cl <- mz2 + 34.9689
    mz2.ch3coo <- mz2 + 59.0133
    mz2_h_h2o <- mz2 - 1.0078 - 18.0106
    mz2.f <- mz2 + 18.9984
    mz2.all.neg <-
      data.frame(mz2_h, mz2.cl, mz2.ch3coo, mz2_h_h2o, mz2.f)
  }

  iden <- rep(NA, length(mz1))
  match.result <- list()

  ## positive
  if (polarity == "positive") {
    for (i in 1:length(mz1)) {
      mz.error <- abs(mz1[i] - mz2.all.pos) * 10 ^ 6 / mz1[i]
      idx <-
        apply(mz.error, 2, function(x) {
          which(x <= mass.tolerance)
        })

      if (length(unlist(idx)) == 0) {
        cat(paste("Feature"), i, "has no matching\n")
        match.result[[i]] <- NA
        next
      }
      else {
        hmdb.name1 <- list()
        hmdb.id1 <- list()
        hmdb.formula1 <- list()
        hmdb.mass1 <- list()
        hmdb.adduct1 <- list()

        add <- c("+H", "+NH4", "+Na", "+H-2H2O", "+H-H2O", "+K")

        for (j in 1:length(idx)) {
          hmdb.name1[[j]] <- HMDB.name[idx[[j]]]
          hmdb.id1[[j]] <- HMDB.ID[idx[[j]]]
          hmdb.formula1[[j]] <- HMDB.formula[idx[[j]]]
          hmdb.mass1[[j]] <- mz2[idx[[j]]]
          hmdb.adduct1[[j]] <- rep(add[j], length(idx[[j]]))
        }

        mass.error1 <-
          lapply(c(1:length(idx)), function(x) {
            mz.error[idx[[x]], x]
          })

        match.result[[i]] <- paste(
          "HMDB.name:",
          unlist(hmdb.name1),
          "HMDB.ID:",
          unlist(hmdb.id1),
          "Formula:",
          unlist(hmdb.formula1),
          "Mass:",
          unlist(hmdb.mass1),
          "HMDB.adduct:" ,
          unlist(hmdb.adduct1),
          "mz.error:",
          unlist(mass.error1)
        )

        rank.order <- order(unlist(mass.error1))
        match.result[[i]] <- match.result[[i]][rank.order]
        if (match.result[[i]] > show) {
          match.result[[i]] <- head(match.result[[i]], show)
        }
        iden[i] <- unlist(hmdb.name1)[rank.order][1]
        match.result[[i]] <- paste(match.result[[i]], collapse = ";")
        cat(paste("Feature"),
            i,
            "has",
            length(rank.order),
            "matching\n")
      }
    }
  }

  ## negative
  if (polarity == "negative") {
    for (i in 1:length(mz1)) {
      mz.error <- abs(mz1[i] - mz2.all.neg) * 10 ^ 6 / mz1[i]
      idx <-
        apply(mz.error, 2, function(x) {
          which(x <= mass.tolerance)
        })

      if (length(unlist(idx)) == 0) {
        cat(paste("Feature"), i, "has no matching\n")
        match.result[[i]] <- NA
        next
      }
      else {
        hmdb.name1 <- list()
        hmdb.id1 <- list()
        hmdb.formula1 <- list()
        hmdb.mass1 <- list()
        hmdb.adduct1 <- list()

        add <- c("-H", "+Cl", "+CH3COO", "-H-H2O", "+F")

        for (j in 1:length(idx)) {
          hmdb.name1[[j]] <- HMDB.name[idx[[j]]]
          hmdb.id1[[j]] <- HMDB.ID[idx[[j]]]
          hmdb.formula1[[j]] <- HMDB.formula[idx[[j]]]
          hmdb.mass1[[j]] <- mz2[idx[[j]]]
          hmdb.adduct1[[j]] <- rep(add[j], length(idx[[j]]))
        }

        mass.error1 <-
          lapply(c(1:length(idx)), function(x) {
            mz.error[idx[[x]], x]
          })

        match.result[[i]] <- paste(
          "HMDB.name:",
          unlist(hmdb.name1),
          "HMDB.ID:",
          unlist(hmdb.id1),
          "Formula:",
          unlist(hmdb.formula1),
          "Mass:",
          unlist(hmdb.mass1),
          "HMDB.adduct:" ,
          unlist(hmdb.adduct1),
          "mz.error:",
          unlist(mass.error1)
        )

        rank.order <- order(unlist(mass.error1))
        match.result[[i]] <- match.result[[i]][rank.order]
        if (match.result[[i]] > show) {
          match.result[[i]] <- head(match.result[[i]], show)
        }
        iden[i] <- unlist(hmdb.name1)[rank.order][1]
        match.result[[i]] <- paste(match.result[[i]], collapse = ";")
        cat(paste("Feature"),
            i,
            "has",
            length(rank.order),
            "matching\n")
      }
    }
  }
  # browser()
  match.result <-
    lapply(match.result, function(x) {
      ifelse(is.null(x), x <- NA, x <- x)
    })
  HMDB.match.result <- unlist(match.result)
  HMDB.identification <- unlist(iden)
  tags <- data.frame(tags, HMDB.match.result, HMDB.identification)
  MetFlowData[["tags"]] <- tags
  MetFlowData[["peak.identification"]] <- "yes"
  return(MetFlowData)
}