
options(warn=2)

library(data.table)
library(stringdist)
library(modeest)

cmdArguments = NULL

main <- function() {
    cmdArguments <<- parseProgramArguments()
    positionalArgs = cmdArguments[[1]]
    do.call(run, as.list(positionalArgs))
}

run <- function(dataFile, reportFile=NULL) {
    if (reportFile == "NA") {
        reportFile = NULL
    }
    stampMapFile = cmdArguments$stampMapFile
    thresholdFile = cmdArguments$thresholdFile
    cat(sprintf("%s Starting analysis ...\n", date()))
    cat(sprintf(" dataFile: %s\n", dataFile))
    cat(sprintf(" reportFile: %s\n", reportFile))
    cat(sprintf(" stampMapFile: %s\n", stampMapFile))
    cat(sprintf(" thresholdFile: %s\n", thresholdFile))
    cat(sprintf(" stampMinReadCount: %s\n", cmdArguments$stampMinReadCount))
    cat(sprintf(" reverseComplementStamps: %s\n", cmdArguments$reverseComplementStamps))
    cat(sprintf(" assignReactionsFromStamps: %s\n", cmdArguments$assignReactionsFromStamps))
    cat(sprintf(" minReadCount: %s\n", cmdArguments$minReadCount))
    cat(sprintf(" minReadPurity: %s\n", cmdArguments$minReadPurity))
    cat(sprintf(" minDistanceThreshold: %s\n", cmdArguments$minDistanceThreshold))
    cat(sprintf(" maxUMILength: %s\n", cmdArguments$maxUMILength))
    cat(sprintf(" outRightThreshold: %s\n", cmdArguments$outRightThreshold))
    if (is.null(cmdArguments$reverseComplementStamps)) {
        cmdArguments$reverseComplementStamps <<- "true"
    }
    if (is.null(cmdArguments$assignReactionsFromStamps)) {
        cmdArguments$assignReactionsFromStamps <<- "false"
    }

    stampFileMap = NULL
    thresholdData = NULL
    thresholdFixed = NA
    if (!is.null(stampMapFile)) {
        cat(sprintf("%s Reading stamp map file %s ...\n", date(), stampMapFile))
        stampFileMap = fread(stampMapFile, header=T, sep="\t")
    }
    if (!is.null(thresholdFile)) {
        cat(sprintf("%s Reading threshold file %s ...\n", date(), thresholdFile))
        thresholdData = fread(thresholdFile, header=T, sep="\t")
    }
    if (!is.null(cmdArguments$minReadCount)) {
        thresholdFixed = as.numeric(cmdArguments$minReadCount)
    }

    cat(sprintf("%s Reading input file %s ...\n", date(), dataFile))
    data = fread(dataFile, header=T, sep="\t")
    cat(sprintf("%s Read %d input records.\n", date(), nrow(data)))

    cat(sprintf("%s Analyzing UMIs ...\n", date()))
    if (!is.null(cmdArguments$minReadPurity)) {
        minPurity = as.numeric(cmdArguments$minReadPurity)
        data = data[data$REPPURITY >= minPurity]
    }
    if (!is.null(cmdArguments$maxUMILength)) {
        maxUMILength = as.numeric(cmdArguments$maxUMILength)
        data = data[nchar(data$RAWUMI) <= maxUMILength]
    }

    totalReadCounts = tapply(data$RAWUMI, data$RAWUMI, length)

    if (asBoolean(cmdArguments$assignReactionsFromStamps) && (0 %in% data$RXN)) {
        stampRxnMap = NULL
        if (is.null(stampFileMap)) {
            cat(sprintf("Warning: assignReactionsFromStamp is set but no stamp files provided\n"))
        } else {
           stampRxnMap = computeStampToRxnMap(stampFileMap)
           if (is.null(stampRxnMap)) {
               cat(sprintf("Warning: Not able to map any stamps unambiguously to reactions\n"))
           }
        }
        if (!is.null(stampRxnMap)) {
            rxn0Data = data[data$RXN == 0]
            mappedRxns = stampRxnMap[rxn0Data$CBC,]$RXN
            mappedRxnIds = stampRxnMap[rxn0Data$CBC,]$RXNID
            mappedRxns[is.na(mappedRxns)] = 0
            data[RXN == 0, RXNID := mappedRxnIds]
            data[RXN == 0, RXN := mappedRxns]
            cat(sprintf("%s Mapped %d out of %d UMIs from rxn 0 using stamps (%1.1f%%)\n",
                        date(), sum(mappedRxns > 0), length(mappedRxns), 100*(sum(mappedRxns > 0)/length(mappedRxns))))
        }
    }

    # Include rxn 0 for QC
    #reactions = sort(unique(data$RXN[!is.na(data$RXN) & data$RXN > 0]))
    reactions = sort(unique(data$RXN[!is.na(data$RXN)]))
    resultTable = NULL

    if (!is.null(cmdArguments$outputReadInfo)) {
        headerData = data[0,]
        write.table(headerData, cmdArguments$outputReadInfo, row.names=F, quote=F, sep="\t", col.names=T)
    }

    for (rxn in reactions) {
        cat(sprintf("%s Analyzing reaction %d ...\n", date(), rxn))
        stamps = NULL
        if (!is.null(stampFileMap)) {
            stampFile = stampFileMap[stampFileMap$RXN == rxn,]$STAMPFILE
            if (length(stampFile) == 0 || is.na(stampFile) || stampFile == "NA") {
                stampFile = NULL
            }
            if (!is.null(stampFile)) {
                stampData = fread(cmd=paste("grep -v ^#",stampFile), header=F)
                cat(sprintf("%s Read stamp data for rxn %d from %s: %d\n", date(), rxn, stampFile, nrow(stampData)))
                stamps = stampData[[1]]
                if (asBoolean(cmdArguments$reverseComplementStamps)) {
                    stamps = reverseComplement(stamps)
                }
            }
        }
        threshold = NA
        umisFromStamps = NULL
        umisFromThreshold = NULL
        rxnData = data[data$RXN == rxn]
        rxnReadCounts = tapply(rxnData$RAWUMI, rxnData$RAWUMI, length)
        consensusRxn = sapply(tapply(data$RXN, data$RAWUMI, consensus), as.integer)
        consensusRxn[is.na(consensusRxn)] = -1
        if (!is.null(stamps)) {
            minReadCount = 0
            if (!is.null(cmdArguments$stampMinReadCount)) {
                minReadCount = as.numeric(cmdArguments$stampMinReadCount)
            }
            umisFromStamps = unique(rxnData$RAWUMI[rxnData$CBC %in% stamps])
            umisFromStamps = umisFromStamps[rxnReadCounts[umisFromStamps] >= minReadCount]
        }
        if (!is.na(thresholdFixed) || !is.null(thresholdData)) {
            threshold = thresholdData$THRESHOLD[thresholdData$RXN == rxn]
            if (is.null(threshold)) {
                threshold = thresholdFixed
            }
            if (!is.na(threshold)) {
                rxnReadCounts = tapply(rxnData$RAWUMI, rxnData$RAWUMI, length)
                umisFromThreshold = names(rxnReadCounts)[rxnReadCounts >= threshold]
            }
        }
        umisFromStamps = umisFromStamps[consensusRxn[umisFromStamps] == rxn]
        umisFromThreshold = umisFromThreshold[consensusRxn[umisFromThreshold] == rxn]
        umis = union(umisFromStamps, umisFromThreshold)
        cat(sprintf("%s Analyzing %d UMIs (%d from STAMPs, %d from read count threshold %s)\n",
                    date(), length(umis), length(umisFromStamps), length(umisFromThreshold), if (is.na(threshold)) "NA" else threshold))
        if (length(umis) == 0) {
            next
        }
        if (length(umis) > 1 && !is.null(cmdArguments$minDistanceThreshold)) {
            minDistThreshold = as.numeric(cmdArguments$minDistanceThreshold)
            # abundance = RXNN
            abundances = tapply(rxnData$RAWUMI, rxnData$RAWUMI, length)[umis]
            filtervec = computeMinDistFilter(umis, abundances, minDistThreshold)
            # df = data.frame(abundances, computeMinDists(umis), filtervec)
            # print(df)
            umis = umis[filtervec]
            cat(sprintf("%s Filtering by mindist %s removed %d UMIs %d -> %d\n",
                        date(), minDistThreshold, sum(!filtervec), length(filtervec), length(umis)))
        }
        if (length(umis) > 0 && !is.null(cmdArguments$outRightThreshold)) {
            outRightThreshold = as.numeric(cmdArguments$outRightThreshold)
            outRight = tapply(rxnData$REPLENGTH, rxnData$RAWUMI, computeOutlierFraction, "R")
            beforeCount = length(umis)
            umis = intersect(umis, names(outRight)[outRight <= outRightThreshold])
            afterCount = length(umis)
            cat(sprintf("%s Filtering by outRight threshold %s removed %d UMIs %d -> %d\n",
                        date(), outRightThreshold, beforeCount-afterCount, beforeCount, afterCount))
        }
        formattedUMIs = formatUMIs(umis)
        umis = umis[order(formattedUMIs)]
        if (length(formattedUMIs) > 0) {
            formattedUMIs = sort(formattedUMIs)
        }
        rxnId = paste(sort(unique(rxnData$RXNID)), collapse=",")
        if (rxnId == "") {
            rxnId = NA
        }
        umiData = rxnData[rxnData$RAWUMI %in% umis]
        umiData$REPUNITS[umiData$REPUNITS <= 0] = NA

        RXN = rep(rxn, length(umis))
        RXNID = rep(rxnId, length(umis))
        RXNN = tapply(rxnData$RAWUMI, rxnData$RAWUMI, length)[umis]
        UMISOURCE = sapply(umis, getUMISource, umisFromStamps, umisFromThreshold)
        NREADS = totalReadCounts[umis]
        READLENGTH = tapply(umiData$READLENGTH, umiData$RAWUMI, median, na.rm=T)[umis]
        CATEGORY = tapply(umiData$CATEGORY, umiData$RAWUMI, consensus)[umis]
        CATEGORYN = tapply(umiData$CATEGORY, umiData$RAWUMI, consensusN)[umis]
        REPLENGTH = tapply(umiData$REPLENGTH, umiData$RAWUMI, halfSampleMode)[umis]
        REPLENGTHN = tapply(umiData$REPLENGTH, umiData$RAWUMI, countNonNA)[umis]
        REPUNITS = tapply(umiData$REPUNITS, umiData$RAWUMI, modeInteger)[umis]
        REPHASCAA = tapply(umiData$REPHASCAA, umiData$RAWUMI, consensus)[umis]
        ALIGNSTART = rep(NA, length(umis))
        ALIGNEND = rep(NA, length(umis))
        ALIGNN = rep(0, length(umis))
        if (all(c("ALIGNSTART","ALIGNEND") %in% names(umiData))) {
            ALIGNN = tapply(umiData$ALIGNSTART, umiData$RAWUMI, countNonNA)[umis]
            ALIGNSTART = suppressWarnings(tapply(umiData$ALIGNSTART, umiData$RAWUMI, modeInteger))[umis]
            ALIGNEND = suppressWarnings(tapply(umiData$ALIGNEND, umiData$RAWUMI, modeInteger))[umis]
        }
        POLYALEN = rep(NA, length(umis))
        if ("POLYALEN" %in% names(umiData)) {
            POLYALEN = tapply(umiData$POLYALEN, umiData$RAWUMI, medianInteger)[umis]
            POLYALEN[CATEGORY != "INTRON1"] = NA
        }
        MINDIST = computeMinDists(umis)[umis]
        OUTLEFT = tapply(umiData$REPLENGTH, umiData$RAWUMI, computeOutlierFraction, "L")
        OUTRIGHT = tapply(umiData$REPLENGTH, umiData$RAWUMI, computeOutlierFraction, "R")

        rxnTable = data.frame(RXN=RXN,
                              RXNID=RXNID,
                              RAWUMI=umis,
                              UMI=formattedUMIs,
                              UMISOURCE=UMISOURCE,
                              NREADS=NREADS,
                              READLENGTH=READLENGTH,
                              RXNN=RXNN,
                              CATEGORY=CATEGORY,
                              CATEGORYN=CATEGORYN,
                              REPLENGTH=REPLENGTH,
                              REPLENGTHN=REPLENGTHN,
                              REPUNITS=REPUNITS,
                              REPHASCAA=REPHASCAA,
                              ALIGNSTART=ALIGNSTART,
                              ALIGNEND=ALIGNEND,
                              ALIGNN=ALIGNN,
                              POLYALEN=POLYALEN,
                              MINDIST=MINDIST,
                              OUTLEFT=OUTLEFT,
                              OUTRIGHT=OUTRIGHT)

        resultTable = rbind(resultTable, rxnTable)
        if (!is.null(cmdArguments$outputReadInfo)) {
            write.table(umiData, cmdArguments$outputReadInfo, row.names=F, quote=F, sep="\t", append=T, col.names=F)
        }
    }
    if (!is.null(reportFile)) {
        reportTable = resultTable
        reportTable$OUTLEFT = sprintf("%1.4f", reportTable$OUTLEFT)
        reportTable$OUTRIGHT = sprintf("%1.4f", reportTable$OUTRIGHT)
        write.table(reportTable, reportFile, row.names=F, quote=F, sep="\t")
    }
    cat(sprintf("%s Analysis done.\n", date()))
}

reverseComplement <- function(seqs) {
    chartr("acgtACGT", "tgcaTGCA", reverseSequences(seqs))
}

reverseSequences <- function(seqs) {
    sapply(lapply(strsplit(seqs, ""), rev), paste, collapse="")
}

getUMISource <- function(umi, umisFromStamps, umisFromThreshold) {
    fromStamps = umi %in% umisFromStamps
    fromThreshold = umi %in% umisFromThreshold
    if (fromStamps && fromThreshold) {
        return("S,T")
    } else if (fromStamps) {
        return("S")
    } else if (fromThreshold) {
        return("T")
    } else {
        return(NA)
    }
}

formatUMIs <- function(umis) {
    sapply(umis, formatUMI)
}

formatUMI <- function(umi) {
    if (nchar(umi) <= 16) {
        return(umi)
    }
    return(sprintf("%s.%s", substring(umi,nchar(umi)-15), substring(umi,1,nchar(umi)-16)))
}

consensus <- function(values) {
    counts = sort(tapply(values,values,length), dec=T)
    if (length(counts) == 0) {
        return(NA)
    }
    if (length(counts) > 1 && counts[[1]] == counts[[2]]) {
        return(NA)
    }
    return(names(counts)[[1]])
}

consensusN <- function(values) {
    counts = sort(tapply(values,values,length), dec=T)
    if (length(counts) == 0) {
        return(0)
    }
    if (length(counts) > 1 && counts[[1]] == counts[[2]]) {
        return(0)
    }
    return(counts[[1]])
}

computeOutlierFraction <- function(values, side) {
    values = values[!is.na(values)]
    if (length(values) == 0) {
        return(NA)
    }
    outliers = countOutliers(values, side=side)
    return(outliers/length(values))
}

countOutliers <- function(replens, mode=NA, side="B", threshold=2.0) {
    if (is.na(mode)) {
        mode = halfSampleMode(replens)
    }
    psd = sqrt(mode)
    if (side == "L") {
        return(sum((replens - mode) < -threshold*psd))
    } else if (side == "R") {
        return(sum((replens - mode) > threshold*psd))
    } else {
        return(sum(abs(replens - mode) > threshold*psd))
    }
}
    
halfSampleMode <- function(values) {
    # Set tie.limit to suppress warnings
    return(hsm(values, tie.limit=1))
}

modeInteger <- function(values) {
    as.integer(round(halfSampleMode(values)))
}

medianInteger <- function(values) {
    as.integer(round(median(values, na.rm=T)))
}

countNonNA <- function(values) {
    return(sum(!is.na(values)))
}

computeMinDists <- function(strings) {
    if (length(strings) == 0) {
        return(NULL)
    }
    if (length(strings) == 1) {
        result = c(NA)
        names(result) = strings
        return(result)
    }
    sdm = stringdistmatrix(strings, strings, method="lv")
    minDists = sapply(1:length(strings), computeMinDist, sdm)
    names(minDists) = strings
    return(minDists)
}

computeMinDist <- function(idx, sdm) {
    cols = setdiff(seq(1,ncol(sdm)), idx)
    v = sdm[idx, cols]
    return(min(v))
}

computeMinDistFilter <- function(umis, abundances, threshold) {
    sdm = stringdistmatrix(umis, umis, method="lv")
    filter = sapply(1:length(umis), computeMinDistFilterHelper, sdm, abundances, threshold)
    return(filter)
}

computeMinDistFilterHelper <- function(idx, sdm, abundances, threshold) {
    cols = setdiff(seq(1,ncol(sdm)), idx)
    v = sdm[idx, cols]
    minDist = min(v)
    if (minDist >= threshold) {
        # cat(sprintf("#DBG: idx %d md=%s t=%s TRUE\n", idx, minDist, threshold))
        return(TRUE)
    }
    ab = abundances[idx]
    mincols = cols[v == minDist]
    maxab = max(abundances[mincols])
    if (ab > maxab || (ab == maxab && idx < min(mincols))) {
        # cat(sprintf("#DBG: idx %d ab=%s mc=%s mab=%s TRUE\n", idx, ab, paste(mincols,collapse=","), max(abundances[mincols])))
        return(TRUE)
    }
    # cat(sprintf("#DBG: idx %d FALSE\n", idx))
    return(FALSE)
}

computeStampToRxnMap <- function(stampFileMap) {
    rxnMap = NULL
    if (length(stampFileMap) > 0 && nrow(stampFileMap) > 0) {
        cat(sprintf("%s Computing reactions from stamps ...\n", date()))
        for (i in 1:nrow(stampFileMap)) {
            rxn = stampFileMap[i,]$RXN
            rxnId = stampFileMap[i,]$RXNID
            stampFile = stampFileMap[i,]$STAMPFILE
            if (length(stampFile) == 0 || is.na(stampFile) || stampFile == "NA") {
                stamps = NULL
                stampFile = NULL
            }
            if (!is.null(stampFile)) {
                stampData = fread(cmd=paste("grep -v ^#",stampFile), header=F)
                cat(sprintf("%s Read stamp data for rxn %d from %s: %d\n", date(), rxn, stampFile, nrow(stampData)))
                stamps = stampData[[1]]
                if (asBoolean(cmdArguments$reverseComplementStamps)) {
                    stamps = reverseComplement(stamps)
                }
            }
            rxnMap = rbind(rxnMap,
                           data.frame(RXN=rep(rxn, length(stamps)),
                                      RXNID=rep(rxnId, length(stamps)),
                                      STAMP=stamps))
        }
        counts = tapply(rxnMap$STAMP, rxnMap$STAMP, length)
        rxnMap = rxnMap[rxnMap$STAMP %in% (names(counts)[counts == 1]),]
        rownames(rxnMap) = rxnMap$STAMP
        cat(sprintf("%s Removed %d/%d non-unique stamps (%1.1f%%), leaving %d [%d] unique stamps\n",
                    date(), sum(counts > 1), length(counts), sum(counts > 1)/length(counts), sum(counts == 1), nrow(rxnMap)))
    }
    return(rxnMap)
}

parseProgramArguments <- function() {
    result = list()
    positional = list()
    result[[""]] = positional
    args = commandArgs()
    if (length(args) == 0) {
        return(result)
    }
    for (i in 1:length(args)) {
        if (i == length(args)) {
            return(result)
        } else if (args[i] == "--args") {
            argpos = i+1
            break
        }
    }
    while (argpos <= length(args)) {
        arg = args[argpos]
        argpos = argpos + 1
        keyword = NULL
        if (nchar(arg) > 2 && substr(arg,1,2) == "--") {
            keyword = substr(arg,3,nchar(arg))
        } else if (nchar(arg) > 1 && substr(arg,1,1) == "-") {
            keyword = substr(arg,2,nchar(arg))
        } else {
            positional = c(positional, arg)
        }
        #cat(sprintf("pos %d kw %s arg %s\n", argpos, keyword, args[argpos]))
        if (!is.null(keyword) && argpos <= length(args)) {
            keyword = as.character(keyword)
            arg = as.character(args[[argpos]])
            argpos = argpos + 1
            result[[keyword]] = c(result[[keyword]], arg)
        }
    }
    result[[1]] = positional
    return(result)
}

asBoolean <- function(arg) {
    if (is.null(arg)) {
        return(FALSE)
    }
    if (is.na(arg)) {
        return(FALSE)
    }
    if (is.logical(arg)) {
        return(arg)
    }
    return(arg == "true")
}

main()
