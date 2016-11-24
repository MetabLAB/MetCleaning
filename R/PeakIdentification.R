#' @title PeakIdentification
#' @description Identify the features using other identification information
#'  data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param path Work directory.
#' @param mz.tolerance: mz tolerance for ms1 and ms2 data matching.
#' @param rt.tolerance: RT tolerance for ms1 and ms2 data matching.
#' @return Return a MetFlowData which has been added peak identification
#' information into tags.
#' @export

### PeakIdentification
PeakIdentification <- function(MetFlowData = MetFlowData,
                               path = "peak identification",
                               identification.information.from = "XCMS",
                               ##parameters for matching
                               peak.number = NULL,
                               mz.tolerance = 30,
                               rt.tolerance = 180,
                               re.match = TRUE) {
  options(warn = -1)
  # browser()
  if (is.null(path)) {
    path <- getwd()
  } else{
    dir.create(path)
  }

  subject <- MetFlowData[["subject"]]
  tags <- MetFlowData[["tags"]]
  if (all(colnames(tags) != "name")) {
    tags[, ncol(tags) + 1] <- paste("feature", c(1:nrow(tags)))
    colnames(tags)[ncol(tags)] <- "name"
  }

  if ("ms2mz" %in% colnames(tags) &
      "ms2rt" %in% colnames(tags)) {
    warning("The data has been done peak identification!")
  }


  if (is.null(peak.number)) {
    peak.number <- nrow(tags)
  }

  # browser()
  peak.number <- as.numeric(peak.number)
  peak.isotopes <- list()
  peak.isotopes[[peak.number + 1]] <- NA
  peak.adduct <- list()
  peak.adduct[[peak.number + 1]] <- NA
  peak.mz <- list()
  peak.mz[[peak.number + 1]] <- NA
  peak.rt <- list()
  peak.rt[[peak.number + 1]] <- NA
  peak.mzerror <- list()
  peak.mzerror[[peak.number + 1]] <- NA
  peak.rterror <- list()
  peak.rterror[[peak.number + 1]] <- NA
  peak.forward <- list()
  peak.forward[[peak.number + 1]] <- NA
  peak.reverse <- list()
  peak.reverse[[peak.number + 1]] <- NA
  peak.name <- list()
  peak.name[[peak.number + 1]] <- NA
  which.file <- list()
  which.file[[peak.number + 1]] <- NA

  text.name <- file.path(path, "Identification.information.txt")
  cat("Identification", file = text.name, append = F)
  # browser()
  ##开始读取二级数据的信息，文件夹中要放入ms2数据，有几个放几个，
  ##命名要写50-300ms2.csv的形式
  data <- dir(path)

  #ms1是一级数据，ms2是二级数据
  ms1 <- tags
  peak <- ms1[, c("mz", "rt", "name")]
  ms2 <- data[grep("ms2", data)]
  ms2.name <- substr(ms2, 1, nchar(ms2) - 4)

  cat(paste("There are ", length(ms2), "ms2 data", "\n"))

  path1 <- file.path(path, "matching result")
  dir.create(path1)

  ##下面开始匹配信息
  if (re.match) {
    cat("Reading ms2 data...")
    ms2 <-
      lapply(ms2, function(x) {
        read.csv(file.path(path, x), stringsAsFactors = FALSE)
      })

    for (i in 1:length(ms2)) {
      msms <- ms2[[i]]
      forward <- as.character(msms[, "hits..forward."])
      reverse <- as.character(msms[, "hits..reverse."])
      forward[is.na(forward)] <- ""
      reverse[is.na(reverse)] <- ""
      #将没有鉴定出来的peak去除
      msms <- msms[forward != '' | reverse != '',]
      cat("\n")
      cat("\n", file = text.name, append = T)
      cat(paste(
        "There are",
        nrow(msms),
        "features are identified in",
        ms2.name[i]
      ))
      cat(
        paste(
          "There are",
          nrow(msms),
          "features are identified in",
          ms2.name[i]
        ),
        file = text.name,
        append = T
      )
      cat("\n", file = text.name, append = T)

      # browser()
      save(msms, file = file.path(path1, ms2.name[i]))
      write.csv(msms, file.path(path1, paste("marker", ms2.name[i], "csv", sep = ".")), row.names = F)
      msmsinfo <- msms[, c("mzmed", "rtmed")]
      #把鉴定出来的peak的信息提取出来
      forward <- as.character(msms[, "hits..forward."])
      reverse <- as.character(msms[, "hits..reverse."])
      forward[is.na(forward)] <- ""
      reverse[is.na(reverse)] <- ""
      name <- as.character(msms[, "name"])
      mz <- as.numeric(msms[, "mzmed"])
      rt <- as.numeric(msms[, "rtmed"])
      adduct <- as.character(msms[, "adduct"])
      isotopes <- as.character(msms[, "isotopes"])
      file <- rep(ms2.name[i], length(forward))

      #开始ms1和ms2数据的匹配
      cat("\n")
      cat(paste("Begin", ms2.name[i], "matching..."))
      cat("\n")
      cat("\n", file = text.name, append = T)

      result <-
        SXTMTmatch(peak,
                   msmsinfo,
                   mz.tolerance = mz.tolerance,
                   rt.tolerance = rt.tolerance)

      if (is.null(result)) {
        next
      }
      else {
        #mz和rt error信息
        mzerror <- result[, "mz error"]
        rterror <- result[, "rt error"]

        #index1是匹配上的ms1的索引，index2是匹配上的ms2的索???
        index1 <- result[, "Index1"]
        index2 <- result[, "Index2"]
        #开始把匹配到的数据信息写入到结果中
        # browser()
        for (i in 1:length(index1)) {
          peak.forward[[index1[i]]] <-
            c(peak.forward[[index1[i]]], forward[index2[i]])
          peak.reverse[[index1[i]]] <-
            c(peak.reverse[[index1[i]]], reverse[index2[i]])
          peak.name[[index1[i]]] <-
            c(peak.name[[index1[i]]], name[index2[i]])
          peak.mz[[index1[i]]] <-
            c(peak.mz[[index1[i]]], mz[index2[i]])
          peak.rt[[index1[i]]] <-
            c(peak.rt[[index1[i]]], rt[index2[i]])
          peak.mzerror[[index1[i]]] <-
            c(peak.mzerror[[index1[i]]], mzerror[i])
          peak.rterror[[index1[i]]] <-
            c(peak.rterror[[index1[i]]], rterror[i])
          peak.isotopes[[index1[i]]] <-
            c(peak.isotopes[[index1[i]]], isotopes[index2[i]])
          peak.adduct[[index1[i]]] <-
            c(peak.adduct[[index1[i]]], adduct[index2[i]])
          which.file[[index1[i]]] <-
            c(which.file[[index1[i]]], file[index2[i]])
        }
      }
    }


    save(
      peak.forward,
      peak.reverse,
      peak.name,
      peak.mz,
      peak.rt,
      peak.mzerror,
      peak.rterror,
      peak.adduct,
      peak.isotopes,
      which.file,
      file = file.path(path1, "msms matching data")
    )
  }

  else {
    load(file.path(path1, "msms matching data"))
  }
  # browser()
  # browser()
  ##这是把所有的MetDDA和LipDDA都弄完之后再开始综
  #下面要对数据进行一个筛选，对于那些一个ms peak对应好几个ms2 peak的情况，根据rterror和mzerror进行判断
  #，只有mzerror的误差在5 ppm之内，才会保留下来，然后在这些peak中选择rterror最小的那个
  #作为该ms1peak应该对应的ms2peak???
  #同时，identification result给每个判断的第一个，也就???
  #打分最高的那个

  path <- getwd()
  path2 <- file.path(path1, "how select one to many")
  dir.create(path2)
  if (all(colnames(ms1) != 'name')) {
    tags[, 1 + ncol(tags)] <- paste("peak", 1:nrow(tags), sep = "")
    colnames(tags)[ncol(tags)] <- "name"
    MetFlowData[["tags"]] <- tags
    ms1 <- tags
    peak <- tags
  }
  ms1name <- as.character(ms1[, "name"])


  for (i in 1:peak.number) {
    if (is.null(peak.name[[i]]) | length(peak.name[[i]]) == 1) {
      next
    }

    else {
      index <-
        which(as.numeric(peak.mzerror[[i]]) - min(as.numeric(peak.mzerror[[i]])) <=
                5)
      index <-
        match(min(as.numeric(peak.rterror[[i]])[index]), as.numeric(peak.rterror[[i]]))
      jpeg(file.path(path2, paste("peak", ms1name[i], ".jpeg", sep = "")))
      plot(
        as.numeric(peak.rterror[[i]]),
        as.numeric(peak.mzerror[[i]]),
        xlab = "rt error",
        ylab = "mz error",
        pch = 20,
        cex.lab = 1.3,
        cex.axis = 1.3,
        xlim = c(0.5 * min(as.numeric(
          peak.rterror[[i]]
        )), 1.5 * max(as.numeric(
          peak.rterror[[i]]
        ))),
        ylim = c(0.5 * min(as.numeric(
          peak.mzerror[[i]]
        )), 1.5 * max(as.numeric(
          peak.mzerror[[i]]
        )))
      )
      points(
        as.numeric(peak.rterror[[i]])[index],
        as.numeric(peak.mzerror[[i]])[index],
        pch = 20,
        col = "red",
        cex = 2
      )
      text(
        as.numeric(peak.rterror[[i]]),
        as.numeric(peak.mzerror[[i]]),
        labels = paste(round(as.numeric(
          peak.rterror[[i]]
        ), 2), round(as.numeric(
          peak.mzerror[[i]]
        ), 2), sep = ","),
        pos = 4
      )
      dev.off()

      peak.forward[[i]] <- peak.forward[[i]][index]
      peak.reverse[[i]] <- peak.reverse[[i]][index]
      peak.name[[i]] <- peak.name[[i]][index]
      peak.mz[[i]] <- peak.mz[[i]][index]
      peak.rt[[i]] <- peak.rt[[i]][index]
      peak.mzerror[[i]] <- peak.mzerror[[i]][index]
      peak.rterror[[i]] <- peak.rterror[[i]][index]
      peak.adduct[[i]] <- peak.adduct[[i]][index]
      peak.isotopes[[i]] <- peak.isotopes[[i]][index]
      which.file[[i]] <- which.file[[i]][index]

    }
  }


  save(
    peak.forward,
    peak.reverse,
    peak.name,
    peak.mz,
    peak,
    rt,
    peak.mzerror,
    peak.rterror,
    peak.adduct,
    peak.isotopes,
    which.file,
    file = file.path(path1, "msms matching data without one to many")
  )


  index <- NULL
  for (i in 1:peak.number) {
    if (is.null(peak.name[[i]])) {
      index <- index
    }
    else {
      index <- c(index, i)
    }
  }


  ms1name <- as.character(peak[, "name"])
  forward <- rep(NA, peak.number)
  reverse <- rep(NA, peak.number)
  ms2name <- rep(NA, peak.number)
  ms2mz <- rep(NA, peak.number)
  ms2rt <- rep(NA, peak.number)
  ms2isotopes <- rep(NA, peak.number)
  ms2adduct <- rep(NA, peak.number)
  mzerror <- rep(NA, peak.number)
  rterror <- rep(NA, peak.number)
  from.file <- rep(NA, peak.number)

  for (i in index) {
    if (is.null(peak.forward[[i]])) {
      forward[i] <- NA
    }
    else {
      forward[i] <- SXTpaste(peak.forward[[i]], sep = "|")
    }

    if (is.null(peak.reverse[[i]])) {
      reverse[i] <- NA
    }
    else {
      reverse[i] <- SXTpaste(peak.reverse[[i]], sep = "|")
    }

    ms2name[i] <- SXTpaste(peak.name[[i]], sep = "|")
    ms2mz[i] <- SXTpaste(peak.mz[[i]], sep = "|")
    ms2isotopes[i] <- SXTpaste(peak.isotopes[[i]], sep = "|")
    ms2adduct[i] <- SXTpaste(peak.adduct[[i]], sep = "|")
    ms2rt[i] <- SXTpaste(peak.rt[[i]], sep = "|")
    mzerror[i] <- SXTpaste(peak.mzerror[[i]], sep = "|")
    rterror[i] <- SXTpaste(peak.rterror[[i]], sep = "|")
    from.file[i] <- SXTpaste(which.file[[i]], sep = "|")

  }


  forward <-
    sapply(forward, function(x) {
      if (is.na(x)) {
        x
      } else {
        ifelse(x == "", NA, x)
      }
    })
  reverse <-
    sapply(reverse, function(x) {
      if (is.na(x)) {
        x
      } else {
        ifelse(x == "", NA, x)
      }
    })
  ms2isotopes <-
    sapply(ms2isotopes, function(x) {
      if (is.na(x)) {
        x
      } else {
        ifelse(x == "", NA, x)
      }
    })
  ms2adduct <-
    sapply(ms2adduct, function(x) {
      if (is.na(x)) {
        x
      } else {
        ifelse(x == "", NA, x)
      }
    })
  names(forward) <-
    names(reverse) <- names(ms2isotopes) <- names(ms2adduct) <- NULL
  # browser()
  lib <- rep(NA, peak.number)
  identification <- rep(NA, peak.number)
  for (i in 1:length(forward)) {
    # cat(i);cat(" ")
    if (!is.na(forward[i]) | !is.na(reverse[i])) {
      if (!is.na(forward[i])) {
        compound <- forward[i]
      }
      else {
        compound <- reverse[i]
      }
      if (regexpr("\\{", compound)[[1]] < 0) {
        lib[i] <- "MetDDA"
        identification[i] <-
          substr(compound, 1, regexpr("Score", compound)[[1]] - 2)
      }
      else {
        lib[i] <- "LipDDA"
        identification[i] <-
          substr(compound,
                 gregexpr("\\{", compound)[[1]][3] + 1,
                 gregexpr("\\}", compound)[[1]][3] - 1)
      }
    }
    else
      (next)
  }
  # browser()
  ###给identification中的重复的内容加上序号
  ide.idx <- which(!is.na(identification))
  ide <- identification[ide.idx]
  dup.ide <- unique(ide[duplicated(ide)])

  if (length(dup.ide) != 0) {
    for (k in 1:length(dup.ide)) {
      temp.idx <- grep(dup.ide[k], ide)
      ide[temp.idx] <-
        paste(dup.ide[k], c(1:length(temp.idx)), sep = "_")
    }

    identification[ide.idx] <- ide
  }
  # browser()
  peak.identification <-
    cbind(
      ms1name,
      peak,
      ms2name,
      ms2mz,
      ms2rt,
      mzerror,
      rterror,
      ms2isotopes,
      ms2adduct,
      forward,
      reverse,
      identification,
      lib,
      from.file
    )
  write.csv(
    peak.identification,
    file.path(path1, "peak.identification.with.many.to.one.csv")
  )



  #下面对于那些多对一的情况进行一下筛选
  marker <-
    peak.identification[!is.na(peak.identification[, "ms2name"]),]
  marker.ms2name <- as.character(marker[, "ms2name"])
  marker.ms2name <- unique(ms2name[!is.na(ms2name)])
  new.marker <- matrix(ncol = ncol(peak.identification) + 1)
  colnames(new.marker) <- c(colnames(marker), "remain")

  path3 <- file.path(path1, "how select many to one")
  dir.create(path3)

  for (i in 1:length(marker.ms2name)) {
    temp <-
      marker[marker[, "ms2name"] == marker.ms2name[i], , drop = FALSE]
    if (nrow(temp) == 1) {
      temp <- cbind(temp, TRUE)
      colnames(temp)[ncol(new.marker)] <- "remain"
      new.marker <- rbind(new.marker, temp)
    }
    else {
      rterror <-
        as.numeric(as.character(temp[, "rterror"]))
      mzerror <- as.numeric(as.character(temp[, "mzerror"]))
      index <- which(mzerror - min(mzerror) <= 5)
      index <- match(min((rterror)[index]), rterror)
      need <- rep(FALSE, nrow(temp))
      need[index] <- TRUE
      temp <- cbind(temp, need)
      colnames(temp)[ncol(new.marker)] <- "remain"
      new.marker <- rbind(new.marker, temp)
      jpeg(file.path(
        path3,
        paste("ms2peak", marker.ms2name[i], "one to many.jpeg")
      ))
      plot(
        rterror,
        mzerror,
        xlab = "rt error",
        ylab = "mz error",
        pch = 20,
        cex.lab = 1.3,
        cex.axis = 1.3,
        xlim = c(0.6 * min(rterror), 1.2 * max(rterror)),
        ylim = c(0.6 * min(mzerror), 1.2 * max(mzerror))
      )
      points(rterror[index],
             mzerror[index],
             pch = 20,
             col = "red",
             cex = 2)
      text(rterror,
           mzerror,
           labels = paste(round(rterror, 2), round(mzerror, 2), sep = ","),
           pos = 4)
      dev.off()
    }
  }

  new.marker <- new.marker[-1, ]
  write.csv(new.marker,
            file.path(path1, "marker.with.many.to.one.csv"))


  ###将那些并不是真正比对上的peak去掉，也就是二级的各个信息改为NA
  # browser()
  remain <- new.marker[, "remain"]
  for (i in 1:length(remain)) {
    if (remain[i]) {
      next
    }
    else {
      new.marker[i, c(5:16)] <- rep(NA, 12)
    }
  }

  write.csv(new.marker,
            file.path(path1, "marker.without.many.to.one.csv"))
  #将鉴定出来的marker信息写入到peak.identification中
  peak.name <- as.character(peak.identification[, "ms1name"])
  marker.name <- as.character(new.marker[, "ms1name"])

  peak.identification[match(marker.name, peak.name), 1:16] <-
    new.marker[, 1:16]
  save(peak.identification, file = file.path(path1, "peak.identification"))
  write.csv(
    peak.identification,
    file.path(path1, "peak.identification.without.many.to.one.csv"),
    row.names = F
  )

  num <- sum(!is.na(peak.identification[,"identification"]))
  cat(
    paste(
      "There are",
      num,
      "features are matched"
    ),
    file = text.name,
    append = T
  )
  cat("\n", file = text.name, append = T)

  # browser()
  if ("ms2mz" %in% colnames(tags) & "ms2rt" %in% colnames(tags)) {
    tags.old <- tags
    tags <- tags[, -c(16:31)]
    tags <- cbind(tags, peak.identification)
    MetFlowData[["tags"]] <- tags
    MetFlowData[["tags.old"]] <- tags.old
    MetFlowData[["peak.identification"]] <- "yes"
  } else {
    tags <- cbind(tags, peak.identification)
    MetFlowData[["tags"]] <- tags
    MetFlowData[["peak.identification"]] <- "yes"
  }
  options(warn = -1)
  return(MetFlowData)
}
