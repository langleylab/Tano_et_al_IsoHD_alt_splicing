library(dplyr)

rankgenomiccoordinates <- function(x) {
    coordinates <- x[,2:3]
    coordinatesvector <- do.call('c', coordinates)
    coordinatesvectorrank <- rank(coordinatesvector, ties.method = 'min')
    coordinatesvectorranknames <- setNames(object = seq_len(length(unique(coordinatesvectorrank))),
                                           nm = sort(unique(coordinatesvectorrank)))
    coordinatesvectorrank <- coordinatesvectorranknames[as.character(coordinatesvectorrank)]
    y <- x
    y$StartCoord <- y$Start
    y$EndCoord <- y$End
    y[,2:3] <- data.frame(matrix(coordinatesvectorrank, nrow = nrow(x)))
    y
}

intersectgenomiccoordinates <- function(ranked) {
    tmp <- list()
    index <- 1

    y <- ranked
    y$intersect <- 0
    for(x in seq_len(nrow(y))) {
        for(i in seq_len(length(tmp))) {
            if((y[x,2]>=tmp[[i]]['start']&&y[x,2]<=tmp[[i]]['end'])||(y[x,3]>=tmp[[i]]['start']&&y[x,3]<=tmp[[i]]['end'])) {
                tmp[[i]]['start'] <- min(tmp[[i]]['start'], y[x, 2])
                tmp[[i]]['end'] <- max(tmp[[i]]['end'], y[x, 3])
                y[x, 'intersect'] <- i
                next   
            }
        }
        if(y[x, 'intersect'] == 0) {
            tmp[[index]] <- c('start' = y[x,2], 'end' = y[x,3])
            y[x, 'intersect'] <- index
            index<-index+1
        }
    }
    y <- split(x = y, f = y$intersect)
    #y <- y[lapply(y,nrow)>1]
    #y <- split(x = y, f = y$intersect)
    y
}

plotgenomiccoordinates <- function(ranked) {
    Strand <- unlist(strsplit(as.character(ranked$Cluster[1]), '_'))[3]
    
    genomiccoordinates.long <- melt(ranked, 
                                    id.vars = c('Start', 'End', 'Junction', 'StartCoord', 'EndCoord'), 
                                    measure.vars = c('Start', 'End'))
    genomiccoordinates.long$Junction <- factor(as.character(genomiccoordinates.long$Junction))
    
    xaxisLabels <- setNames(object = unique(sort(genomiccoordinates.long$value)),
                            nm = unique(sort(do.call('c', genomiccoordinates.long[,c('StartCoord', 'EndCoord')]))))
    
    exons <- rbind(data.frame(Junction = rep(genomiccoordinates.long$Junction[genomiccoordinates.long$variable == 'Start'], times = 2), 
                              value = c(genomiccoordinates.long$value[genomiccoordinates.long$variable == 'Start'], genomiccoordinates.long$value[genomiccoordinates.long$variable == 'Start'] - .5),
                              group = 'Start'),
                   data.frame(Junction = rep(genomiccoordinates.long$Junction[genomiccoordinates.long$variable == 'End'], times = 2), 
                              value = c(genomiccoordinates.long$value[genomiccoordinates.long$variable == 'End'], genomiccoordinates.long$value[genomiccoordinates.long$variable == 'End'] + .5),
                              group = 'End'))
    exonLabel <- data.frame(Junction = rep(levels(exons$Junction), times = 2),
                            value = c(tapply(X = exons$value[exons$group == 'Start'], INDEX = exons$Junction[exons$group == 'Start'], FUN = mean),
                                      tapply(X = exons$value[exons$group == 'End'], INDEX = exons$Junction[exons$group == 'End'], FUN = mean)),
                            label = ifelse(Strand == '+', '>', '<'))
    
    ggplot(data = genomiccoordinates.long,
           mapping = aes(x = value, y = Junction)) +
    theme_bw() +
    geom_line(mapping = aes(x = value, y = Junction), linetype = 'dashed') +
    scale_x_continuous(name = '', breaks = xaxisLabels, labels = names(xaxisLabels)) +
    geom_line(mapping = aes(x = value, y = Junction, group = paste0(Junction,':',group)), data = exons, size = 5) +
    geom_text(mapping = aes(x = value, y = Junction, label = label, color = 'white'), data = exonLabel, size = 5, nudge_y = 0.05) +
    scale_color_manual(values = c('white' = 'white'), guide = F) +
    labs(title = paste0(ranked[1,'Cluster'], '.', ranked[1,'intersect']))
}

getmatchingjunctions <- function(a,b) {
    if(nrow(a)==0) return(FALSE)
    if(nrow(b)==0) return(FALSE)
    unlist(apply(a,
                 1,
                 function(x) {
                     any(unlist(apply(b,
                                      1,
                                      function(y) {
                                          as.numeric(x[1]) == as.numeric(y[1]) &
                                          as.numeric(x[2]) == as.numeric(y[2])
                                      })))
                 }))
}

plotgenomiccoordinateswithexpression <- function (ranked,plot.title=NULL) 
{
    Strand <- unlist(strsplit(as.character(ranked$Cluster[1]), 
        "_"))[3]

    chromsome <- tryCatch(expr={ranked$Chr[1]},error=function(cond){return(NULL)})
    genes <- tryCatch(expr={ranked$genes[1]},error=function(cond){return(NULL)})
    
    genomiccoordinates.long <- melt(
        ranked, 
        id.vars = c("Start", "End","Junction","StartCoord","EndCoord",'deltapsi','padj'), 
        measure.vars = c("Start", "End")
    ) %>% 
    data.frame
    
    genomiccoordinates.long$Junction <- factor(as.character(genomiccoordinates.long$Junction))
    xaxisLabels <- setNames(object = unique(sort(genomiccoordinates.long$value)), 
        nm = unique(sort(do.call("c", genomiccoordinates.long[, 
            c("StartCoord", "EndCoord")]))))
    exons <- rbind(data.frame(Junction = rep(genomiccoordinates.long$Junction[genomiccoordinates.long$variable == 
        "Start"], times = 2), value = c(genomiccoordinates.long$value[genomiccoordinates.long$variable == 
        "Start"], genomiccoordinates.long$value[genomiccoordinates.long$variable == 
        "Start"] - 0.5), group = "Start"), data.frame(Junction = rep(genomiccoordinates.long$Junction[genomiccoordinates.long$variable == 
        "End"], times = 2), value = c(genomiccoordinates.long$value[genomiccoordinates.long$variable == 
        "End"], genomiccoordinates.long$value[genomiccoordinates.long$variable == 
        "End"] + 0.5), group = "End"))
    exonLabel <- data.frame(Junction = rep(levels(exons$Junction), 
        times = 2), value = c(tapply(X = exons$value[exons$group == 
        "Start"], INDEX = exons$Junction[exons$group == "Start"], 
        FUN = mean), tapply(X = exons$value[exons$group == "End"], 
        INDEX = exons$Junction[exons$group == "End"], FUN = mean)), 
        label = ifelse(Strand == "+", ">", "<"))

    genomiccoordinates.long$Junction <- genomiccoordinates.long$Junction %>% factor
    #return(genomiccoordinates.long)
    g <- ggplot(
        data = genomiccoordinates.long, 
        mapping = aes(
            x = value, 
            y = Junction,
            color=ifelse(deltapsi>0,'up','down') %>% factor(levels=c('up','down'),labels=c('up-regulated','down-regulated'))
        )
    ) + 
    geom_line(
        mapping = aes(
            x = value, 
            y = Junction
        ), 
        size=.5,linetype = "dashed"
    ) + 
    scale_x_continuous(
        name=chromsome, 
        breaks = xaxisLabels, labels = names(xaxisLabels)
    ) + 
    geom_line(
        mapping=aes(
            x=value, 
            y=Junction, 
            group=paste0(Junction,":", group)), 
        data=exons,size=5,color='black'
    ) + 
    geom_text(
        mapping = aes(
            x = value, 
            y = Junction, 
            label = label
        ), 
        data=exonLabel, 
        size=5,nudge_y=0,color="white"
    ) + 
    scale_color_manual(
        name='ΔPSI',
        values=c('up-regulated'='coral','down-regulated'='forestgreen')
    ) + 
    labs(title = paste0(ranked[1, "Cluster"],".",ranked[1, "intersect"],' [',genes,'] - ',plot.title)) +
    theme_bw()
    
    if (ranked$ScaledExp %>% is.list) {
        ranked$ScaledExp <- ranked$ScaledExp %>% setNames(nm=ranked$Junction)
        
        expdata <- ranked$ScaledExp %>% names %>% 
        lapply(function(Junction) {
            ranked$ScaledExp[[Junction]] %>% names %>% 
            lapply(function(Condition) {
                data.frame(
                    Junction=Junction,
                    Condition=Condition,
                    Expression=ranked$ScaledExp[[Junction]][[Condition]]
                )
            }) %>% 
            (function(dflist) do.call(rbind,dflist))
        }) %>% 
        (function(dflist) do.call(rbind,dflist))
        
        expdata$Junction <- expdata$Junction %>% factor(levels=genomiccoordinates.long$Junction %>% levels)
        expdata$Expression <- expdata$Expression + (expdata$Junction %>% as.numeric)
        expdata$value <- (tapply(genomiccoordinates.long$value,genomiccoordinates.long$Junction,mean))[expdata$Junction]
    
        #return(expdata)
        g <- g +
        geom_jitter(
            data=expdata,
            mapping=aes(
                x=value,
                y=Expression,
                fill=Condition
            ),
            color='black',shape=21,stroke=0,
            height=0,width=.2
        )
    }
    g
}

# split non-overlapping clusters
renameClusterLeafcutterIntronsList <- function(ds_effect_sizes=NULL) {
    if (ds_effect_sizes %>% is.null) return(NULL)
    
    junctionlist <- effect_sizes %>% split(f=ds_effect_sizes$Cluster)
    
    junctionintersectlist <- junctionlist %>% lapply(function(x) x %>% rankgenomiccoordinates %>% intersectgenomiccoordinates)
    
    junctionintersectlong <- do.call('c',junctionintersectlist)
    
    junctionintersectrenamed <- junctionintersectrenamed %>% lapply(function(x) {
        y <- x[,2:3] %>% unlist
        rename <- seq_len(length(unique(sort(y)))) %>% 
        setNames(nm=unique(sort(y)))
        x[,2:3] <- matrix(data = rename[as.character(y)], nrow = nrow(x))
        return(x)
    })
    
    return(junctionintersectrenamed)
}

# preliminary annotation of putative AS event type
annotatePutativeASEventTypeList <- function(event_list=NULL) {
    if (event_list %>% is.null) return(NULL)
    
    junctionASevents <- rep('',times=event_list %>% length) %>% 
        setNames(nm=event_list %>% names)
    
    # lone junctions
    junctionASevents[(junctionASevents %>% lapply(nrow))==1] <- 'intron'
    
    # Alt 5' ss (A5SS) forward genes
    A5SSformatrix <- matrix(c(1,2,3,3), nrow = 2)
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(dim(y[,2:3]) == dim(A5SSformatrix) & all(y[,2:3] == A5SSformatrix))
        }))==TRUE & grepl(pattern = '\\+', names(x))
        ] <- 'A5SS'

    # Alt 5' ss (A5SS) reverse genes
    A5SSrevmatrix <- matrix(c(1,1,2,3), nrow = 2)
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(dim(y[,2:3]) == dim(A5SSrevmatrix) & all(y[,2:3] == A5SSrevmatrix))
        }))==TRUE & grepl(pattern = '\\-', names(x))
        ] <- 'A5SS'
    
    # Alt 3' ss (A3SS) forward genes
    A3SSformatrix <- matrix(c(1,1,2,3), nrow = 2)
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(dim(y[,2:3]) == dim(A3SSformatrix) & all(y[,2:3] == A3SSformatrix))
        }))==TRUE & grepl(pattern = '\\+', names(x))
        ] <- 'A3SS'
    
    # Alt 3' ss (A3SS) reverse genes
    A3SSrevmatrix <- matrix(c(1,2,3,3), nrow = 2)
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(dim(y[,2:3]) == dim(A3SSrevmatrix) & all(y[,2:3] == A3SSrevmatrix))
        }))==TRUE & grepl(pattern = '\\-', names(x))
        ] <- 'A3SS'
    
    # Mixed alt ss
    MSSmatrix <- matrix(c(1,2,4,3), nrow = 2)
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(dim(y[,2:3]) == dim(MSSmatrix) & all(y[,2:3] == MSSmatrix))
        }))==TRUE & grepl(pattern = '[\\+\\-]', names(x))
        ] <- 'MSS'
    
    # Skipped exons
    SEmatrix <- matrix(c(1,1,3,2,4,4), nrow = 3)
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(dim(y[,2:3]) == dim(SEmatrix) & all(y[,2:3] == SEmatrix))
        }))==TRUE & grepl(pattern = '[\\+\\-]', names(x))
        ] <- 'SE'
    
    # multi 5' ss forward
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(nrow(y)>2 & length(unique(y[,2]))==1)
        }))==TRUE & grepl(pattern = '\\+', names(x))
        ] <- 'A5SSmulti'
    
    # multi 5' ss reverse
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(nrow(y)>2 & length(unique(y[,3]))==1)
        }))==TRUE & grepl(pattern = '\\-', names(x))
        ] <- 'A5SSmulti'
    
    # multi 3' ss reverse
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(nrow(y)>2 & length(unique(y[,2]))==1)
        }))==TRUE & grepl(pattern = '\\-', names(x))
        ] <- 'A3SSmulti'

    # multi 3' ss forward
    junctionASevents[
        (junctionASevents %>% lapply(function(y) {
            return(nrow(y)>2 & length(unique(y[,3]))==1)
        }))==TRUE & grepl(pattern = '\\+', names(x))
        ] <- 'A3SSmulti'
}

annotatePutativeASEventTypeTable <- function(event_list=NULL) {
    if (event_list %>% is.null) return(NULL)
    
    intronjunctionsAS_list <- event_list %>% annotatePutativeASEventTypeList
    
    intronjunctionsAS_table <- data.frame(matrix(nrow=event_list %>% length,ncol=0),row.names=x %>% names) %>% 
        mutate(NumberOfJunctions=event_list %>% lapply(nrow) %>% unlist)
    
    # Skipped exons
    # pair of junctions with one "end" greater than another "start" by exactly one in terms of ranked position
    clusternumskippedexons <- event_list %>% 
        lapply(function(cluster) {
            nrow(unique(do.call(rbind,apply(
                cluster,
                1,
                function(x) {
                    df <- data.frame(matrix(nrow = 0, ncol = 2))
                    df.list <- apply(cluster,
                                     1,
                                     function(y) {
                                         if(as.numeric(x[2])-1==as.numeric(y[3])) return(data.frame(matrix(c(as.numeric(y[3]),as.numeric(x[2])), nrow = 1)))
                                         else return(data.frame(matrix(nrow = 0, ncol = 2)))
                                     })
                    do.call(rbind,df.list)
                }
            ))))
        }) %>% 
        unlist
    intronjunctionsAS_table <- intronjunctionsAS_table %>% 
        mutate(
            NumberOfSkippedExons = clusternumskippedexons,
            SkippedExons = clusternumskippedexons>0
        )
    
    # Cassette exons
    # containing at least one pair of junctions with one "end" greater than another "start" by exactly one in terms of ranked position, and exactly covered by exon spanning reads skipping the whole cassette
    clustersnumcassetteexons <- event_list %>% 
        lapply(function(cluster) {
            cassettejunctions <- unique(do.call(rbind,apply(
                cluster,
                1,
                function(x) {
                    df <- data.frame(matrix(nrow = 0, ncol = 2))
                    df.list <- apply(cluster,
                                     1,
                                     function(y) {
                                         if(as.numeric(x[2])-1==as.numeric(y[3])) return(data.frame(matrix(c(as.numeric(y[2]),as.numeric(x[3])), nrow = 1)))
                                         else return(data.frame(matrix(nrow = 0, ncol = 2)))
                                     })
                    do.call(rbind,df.list)
                }
            )))
            sum(getmatchingjunctions(cassettejunctions, cluster[,2:3]))
        }) %>% 
        unlist
    intronjunctionsAS_table <- intronjunctionsAS_table %>% 
        mutate(
            NumberOfCassetteExons = clustersnumcassetteexons,
            CassetteExons = clustersnumcassetteexons>0
        )

    # Alternative 5' splice site (A5SS)
    # containing at least one pair of junctions with "start" the first ('1') and second ('2') positions in '+' clusters, and "end" last (max) and second (max-1) last positions in '-' clusters
    # A5SSFor <- event_list[
    #     event_list %>% 
    #         lapply(function(x) {
    #             any(x[,2]==1) & 
    #                 any(x[,2]==2) & 
    #                 nrow(x)==2})==T &
    #         grepl(pattern = '\\+', x = event_list %>% names)
    #     ]
    A5SSFor <- event_list[
        intronjunctionsAS_list == 'A5SS' & grepl(pattern = '\\+', x = names(intronjunctionsAS_list))
    ]
    # A5SSRev <- event_list[
    #     event_list %>% 
    #         lapply(function(x) {
    #             any(x[,3]==(max(x[,3])-1)) && 
    #                 any(x[,3]==max(x[,3])) && 
    #                 nrow(x)==2})==T &
    #         grepl(pattern = '\\-', x = event_list %>% names)
    #     ]
    A5SSRev <- event_list[
        intronjunctionsAS_list == 'A5SS' & grepl(pattern = '\\-', x = names(intronjunctionsAS_list))
    ]
    intronjunctionsAS_table <- intronjunctionsAS_table %>% 
        mutate(
            A5SS = FALSE
        )
    intronjunctionsAS_table$A5SS[
        event_list %>% lapply(function(x) {
            any(x[,2]==2) & any(x[,2]==1) 
        }) & grepl(pattern = '\\+', x = event_list %>% names)
    ] <- TRUE
    intronjunctionsAS_table$A5SS[
        event_list %>% lapply(function(x) {
            any(x[,3]==(max(x[,3])-1)) & any(x[,3]==max(x[,3]))
        }) & grepl(pattern = '\\-', x = event_list %>% names)
    ] <- TRUE    
    
    # Alternative 3' splice site (A3SS)
    # containing at least one pair of junctions with "start" the first ('1') and second ('2') positions in '-' clusters, and "end" last (max) and second (max-1) last positions in '+' clusters
    # A3SSRev <- event_list[
    #     event_list %>% 
    #         lapply(function(x) {
    #             any(x[,2]==1) & 
    #               any(x[,2]==2) & 
    #               nrow(x)==2})==T &
    #         grepl(pattern = '\\+', x = event_list %>% names)
    #     ]
    A3SSRev <- event_list[
        intronjunctionsAS_list == 'A3SS' & grepl(pattern = '\\-', x = names(intronjunctionsAS_list))
    ]
    # A3SSFor <- event_list[
    #     event_list %>% 
    #         lapply(function(x) {
    #             any(x[,3]==(max(x[,3])-1)) & 
    #                 any(x[,3]==max(x[,3])) & 
    #                 nrow(x)==2})==T &
    #         grepl(pattern = '\\+', x = event_list %>% names)
    #     ]
    A3SSFor <- event_list[
        intronjunctionsAS_list == 'A3SS' & grepl(pattern = '\\+', x = names(intronjunctionsAS_list))
    ]
    intronjunctionsAS_table <- intronjunctionsAS_table %>% 
        mutate(
            A3SS = FALSE
        )
    intronjunctionsAS_table$A3SS[
        event_list %>% lapply(function(x) {
            any(x[,2]==2) & any(x[,2]==1) 
        }) & grepl(pattern = '\\-', x = event_list %>% names)
    ] <- TRUE
    intronjunctionsAS_table$A3SS[
        event_list %>% lapply(function(x) {
            any(x[,3]==(max(x[,3])-1)) & any(x[,3]==max(x[,3]))
        }) & grepl(pattern = '\\+', x = event_list %>% names)
    ] <- TRUE
    
    # Multiple skipped exons
    # contains junctions skipping multiple skipped exons and are observed
    clustersnummultiskippedexons <- event_list %>% 
        lapply(function(cluster) {
            skipped_exons <- data.frame(matrix(nrow = 0, ncol = 4))
            
            for(x in seq_len(nrow(cluster))) {
                for(y in seq_len(nrow(cluster))) {
                    if(as.numeric(cluster[x,2])-1 == as.numeric(cluster[y,3])) {
                        skipped_exons <- rbind(skipped_exons,c(cluster[y,3],cluster[x,2],cluster[y,2],cluster[x,3]))
                    }
                }
            }
            
            if(nrow(skipped_exons)==0) return(0)
            
            skipped_exons <- skipped_exons %>% unique
            colnames(skipped_exons) <- c('Start','End','BoundStart','BoundEnd')
            
            expanded.skipped_exons <- expand.grid(
                BoundStart = skipped_exons$BoundStart,
                BoundEnd = skipped_exons$BoundEnd
            )
            expanded.skipped_exons <- expanded.skipped_exons[
                expanded.skipped_exons$BoundEnd>(expanded.skipped_exons$BoundStart+4),
            ]
            expanded.skipped_exons <- expanded.skipped_exons %>% unique
            obs.skipped_exons <- skipped_exons[
                getmatchingjunctions(skipped_exons[,3:4],cluster[,2:3]),
            ]
            obs.expanded.skipped_exons <- expanded.skipped_exons[
                getmatchingjunctions(expanded.skipped_exons,cluster[,2:3]),
            ]
            obs.expanded.skipped_exons$numofskippedexons <- obs.expanded.skipped_exons %>% 
                apply(1,function(x) {
                    sum(unique(skipped_exons[,1:2])$Start > x[1] &
                            unique(skipped_exons[,1:2])$End < x[2])
                })
            obs.expanded.skipped_exons$numofcassetteexons <- obs.expanded.skipped_exons %>% 
                apply(1,function(x) {
                    sum(obs.skipped_exons[,3] == x[1] &
                            obs.skipped_exons[,4] == x[2])
                })
            sum(obs.expanded.skipped_exons$numofcassetteexons < obs.expanded.skipped_exons$numofskippedexons)
        }) %>% unlist
    
    clustersmultiskippedexonsnotcompletelycoveredchained <- event_list %>% 
        lapply(function(cluster) {
            cluster$Junction <- as.character(cluster$Junction)
            
            # find all skipped exons and return the junctions and boundaries
            skipped_exons <- data.frame(matrix(nrow = 0, ncol = 8))
            for(x in seq_len(nrow(cluster))) {
                for(y in seq_len(nrow(cluster))) {
                    if(as.numeric(cluster[x,2])-1 == as.numeric(cluster[y,3])) {
                        tmp.df <- data.frame(cluster[y,3],cluster[x,2],cluster[y,2],cluster[x,3],
                                             cluster[y,11],cluster[x,10])
                        tmp.df$JunctionStart <- cluster[y,4]
                        tmp.df$JunctionEnd <- cluster[x,4]
                        skipped_exons <- rbind(skipped_exons,tmp.df)
                    }
                }
            }
            if(nrow(skipped_exons)==0) return(F)
            skipped_exons <- unique(skipped_exons)
            colnames(skipped_exons) <- c('Start','End','BoundStart','BoundEnd','CoordStart','CoordEnd','JunctionStart','JunctionEnd')
            
            # predict possible multiple skipped exons junctions
            expanded.skipped_exons <- expand.grid(BoundStart = skipped_exons$BoundStart,
                                                  BoundEnd = skipped_exons$BoundEnd)
            expanded.skipped_exons <- expanded.skipped_exons[expanded.skipped_exons$BoundEnd>(expanded.skipped_exons$BoundStart+4),]
            expanded.skipped_exons <- unique(expanded.skipped_exons)
            
            # check if skipped exons are cassettee exons
            obs.skipped_exons <- skipped_exons[getmatchingjunctions(skipped_exons[,3:4],cluster[,2:3]),]
            
            # match predicted multiple-skipped-exon junctions with experimentally observed junctions
            obs.expanded.skipped_exons <- expanded.skipped_exons[getmatchingjunctions(expanded.skipped_exons,cluster[,2:3]),]
            obs.expanded.skipped_exons$matchedjunction <- apply(obs.expanded.skipped_exons,1,function(foo) {
                cluster$Junction[cluster[2]==foo[1]&cluster[3]==foo[2]][1]
            })
            
            # retain only chained exons
            any(unlist(lapply(1:nrow(obs.expanded.skipped_exons),function(x) {
                y <- obs.expanded.skipped_exons[x,]
                #skipped_exons.contained <- skipped_exons[skipped_exons$BoundStart>=y[1] & skipped_exons$BoundEnd<=y[2],]
                skipped_exons$index <- 1:nrow(skipped_exons)
                
                forwardchains <- c(y[1])
                forwardchainexons <- c()
                for(i in order(skipped_exons$BoundStart, decreasing = F)) {
                    if(skipped_exons[i,3] %in% forwardchains) {
                        forwardchains <- c(forwardchains, skipped_exons[i,2])
                        forwardchainexons <- c(forwardchainexons, skipped_exons[i,9])
                    }
                }
                
                reversechains <- c(y[2])
                reversechainexons <- c()
                for(i in order(skipped_exons$BoundEnd, decreasing = T)) {
                    if(skipped_exons[i,4] %in% reversechains) {
                        reversechains <- c(reversechains, skipped_exons[i,1])
                        reversechainexons <- c(reversechainexons, skipped_exons[i,9])
                    }
                }
                
                skipped_exons.chained <- skipped_exons[skipped_exons$index %in% forwardchainexons & skipped_exons$index %in% reversechainexons,]
                skipped_exons.chained.notcassette <- skipped_exons.chained[!getmatchingjunctions(skipped_exons.chained[,1:2],obs.skipped_exons[,1:2]),]
                nrow(skipped_exons.chained.notcassette)
            })>1))
        }) %>% unlist
    intronjunctionsAS_table <- intronjunctionsAS_table %>% 
        mutate(
            NumberOfMultiSkippedExons = clustersnummultiskippedexons,
            MultiSkippedExons = clustersnummultiskippedexons>0
        )
    
    # Mutually exclusive exons
    clustersnummutuallyexclusiveexons <- event_list %>% 
        lapply(function(cluster) {
            uniqueexons <- unique(
                do.call(
                    rbind,cluster %>%
                        apply(1,function(x) {
                            df <- data.frame(matrix(nrow = 0, ncol = 2))
                            df.list <- apply(cluster,
                                             1,
                                             function(y) {
                                                 if(as.numeric(x[2])-1==as.numeric(y[3])) return(data.frame(matrix(c(as.numeric(y[3]),as.numeric(x[2])), nrow = 1)))
                                                 else return(data.frame(matrix(nrow = 0, ncol = 2)))
                                             })
                            do.call(rbind,df.list)
                        })))
            #if (nrow(uniqueexons)>=2) {
            betweenexonjunctions <- expand.grid(Start = uniqueexons[,2], End = uniqueexons[,1])
            betweenexonjunctions <- betweenexonjunctions[order(betweenexonjunctions$Start,betweenexonjunctions$End),]
            betweenexonjunctions <- subset(betweenexonjunctions, Start<=End)
            betweenexonjunctions <- betweenexonjunctions[!duplicated(betweenexonjunctions$Start),]
            if(nrow(betweenexonjunctions)>0) {
                sum(!getmatchingjunctions(betweenexonjunctions, cluster[,2:3]))
            } else {
                return(F)
            }              
            #} else {
            #    return(F)
            #}
        }) %>% unlist
    intronjunctionsAS_table <- intronjunctionsAS_table %>% 
        mutate(
            NumberOfMutuallyExclusiveExonPairs = clustersnummutuallyexclusiveexons,
            MutuallyExclusiveExons = NumberOfMutuallyExclusiveExonPairs>0
        )
    
    intronjunctionsAS_table
}

annotatePutativeASEventType <- function(eventtype_table=NULL) {
    if (eventtype_table %>% is.null) return(NULL)
    
    eventtype_table %>% 
        mutate(
            ASType=apply(1,function(x) {
                if ((c(x[['CassetteExons']],x[['A5SS']],x[['A3SS']],x[['MultiSkippedExons']],x[['MutuallyExclusiveExons']]) %>% sum)>=2) return ('Mixed')
                else if (x[['NumberOfJunctions']]==3 & x[['NumberOfSkippedExons']]==1 & x[['SkippedExons']] & x[['NumberOfCassetteExons']]==1 & x[['CassetteExons']] & !x[['A5SS']] & !x[['A3SS']] & !x[['MultiSkippedExons']] & !x[['MutuallyExclusiveExons']]) return ('CassetteExon')
                else if (x[['NumberOfJunctions']]==2 & !x[['SkippedExons']] & !x[['CassetteExons']] & x[['A5SS']] & !x[['A3SS']] & !x[['MultiSkippedExons']] & !x[['MutuallyExclusiveExons']]) return ('A5SS')
                else if (x[['NumberOfJunctions']]==2 & !x[['SkippedExons']] & !x[['CassetteExons']] & !x[['A5SS']] & x[['A3SS']] & !x[['MultiSkippedExons']] & !x[['MutuallyExclusiveExons']]) return ('A3SS')
                else if (x[['NumberOfJunctions']]==4 & !x[['CassetteExons']] & !x[['A5SS']] & !x[['A3SS']] & !x[['MultiSkippedExons']] & x[['MutuallyExclusiveExons']]) return ('MXE')
                else if (x[['NumberOfJunctions']]==1) return ('PossibleRI')
                else return ('Unknown')
            })
        )
}

# match cassette exons
.match_cassetteexon <- function(
    intronjunctionsdf=NULL,ASType='CassetteExon',
    EVENT_INFO=EVENT_INFO.hg38.tab,EVENT_CONSERVATION=EVENT_CONSERVATION.tab.hg38mm10,
    deltapsicutoff=.01,padjcutoff=.1,containnoveljunction=FALSE,alt=FALSE
    ) {
        #message(deltapsicutoff)
        novelcassetteexon <- intronjunctionsdf %>% 
        subset(ClusterASType==ASType) %>% 
        (function(mydf) {
            #message(deltapsicutoff)
            mydflist <- mydf %>% 
            split(f=mydf$renamedCluster)

            mydflistinclude <- mydflist %>% 
            lapply(
                function(mydf,deltapsicutoff.=deltapsicutoff,containnoveljunction.=containnoveljunction) {
                    mydf %>% 
                    with((((deltapsi %>% abs)>=deltapsicutoff.) %>% any)&(clusterContainsNovelJunction %>% any)==containnoveljunction.)
                }
            ) %>% 
            unlist
            if ((mydflistinclude %>% sum)==0) {
                return(NULL)
            }

            outdf <- mydflist[mydflistinclude] %>% 
            lapply(function(x) {
                coordfull <- paste0(x %>% subset(Start==1&End==4) %>% with(StartCoord),',',x %>% subset(Start==1&End==2) %>% with(EndCoord),'-',x %>% subset(Start==3&End==4) %>% with(StartCoord),',',x %>% subset(Start==1&End==4) %>% with(EndCoord))
                coordstartend <- paste0(x %>% subset(Start==1&End==2) %>% with(EndCoord),'-',x %>% subset(Start==3&End==4) %>% with(StartCoord))
                x$coordfull <- gsub(' *','',coordfull)
                x$coordstartend <- paste0(x$Chr,':',gsub(' *','',coordstartend))
                x
            }) %>% 
            (function(dflist) do.call(rbind,dflist))

            if (alt) {
                outdf2 <- outdf %>% 
                merge(
                    y=EVENT_INFO %>% subset(COMPLEX=='S'|COMPLEX=='S*'|COMPLEX=='C1'|COMPLEX=='C1*'|COMPLEX=='A_S'|COMPLEX=='A_C1'),
                    #by.x='genes',by.y='GENE',all.x=TRUE
                    by.x='coordstartend',by.y='COORD_o',all.x=TRUE
                )
            } else {
                outdf2 <- outdf %>% 
                merge(
                    y=EVENT_INFO %>% subset(COMPLEX=='S'|COMPLEX=='S*'|COMPLEX=='C1'|COMPLEX=='C1*'),
                    #by.x='genes',by.y='GENE',all.x=TRUE
                    by.x='coordstartend',by.y='COORD_o',all.x=TRUE
                )
            }
           
            # message(outdf$EVENT %>% is.na %>% sum)
            outdf2
        })

        if (novelcassetteexon %>% is.null) {
            message(paste0('0 cassette exons'))
            return(NULL)
        }
    
        outdf <- novelcassetteexon %>% 
        distinct(Chr,JuncStart=StartCoord,JuncEnd=EndCoord,Junction,renamedCluster,genes,splice_site,anchor,known_junction=(known_junction=='1'),ClusterASType,deltapsi,padj,clusterContainsNovelJunction,EVENT,COORD_o=coordstartend,LE_o=gsub(' *','',LE_o) %>% as.numeric,FULL_CO,COMPLEX,REF_CO,LE_n=gsub(' *','',LE_n) %>% as.numeric,CO_C1,CO_A,CO_C2,Seq_C1,Seq_A,Seq_C2) %>% 
        (function(mydf) mydf[,c('Chr','JuncStart','JuncEnd','Junction','renamedCluster','genes','splice_site','anchor','known_junction','ClusterASType','deltapsi','padj','clusterContainsNovelJunction','EVENT','COORD_o','LE_o','FULL_CO','COMPLEX','REF_CO','LE_n','CO_C1','CO_A','CO_C2','Seq_C1','Seq_A','Seq_C2')]) %>% 
        subset((padj %>% as.numeric)<=padjcutoff) %>% 
        merge(EVENT_CONSERVATION,by.x='EVENT',by.y='EventID',all.x=TRUE)

        message(paste0(outdf %>% with(renamedCluster) %>% unique %>% length,' cassette exons'))
        return(outdf)
}

# match alt. 5' splice site
.match_alt5ss <- function(
    intronjunctionsdf=NULL,ASType='A5SS',
    EVENT_INFO=EVENT_INFO.hg38.tab,EVENT_CONSERVATION=EVENT_CONSERVATION.tab.hg38mm10,
    deltapsicutoff=.01,padjcutoff=.1,containnoveljunction=FALSE
    ) {
        novela5ss <- intronjunctionsdf %>% 
        subset(ClusterASType==ASType) %>% 
        (function(mydf) {
            mydflist <- mydf %>% 
            split(f=mydf$renamedCluster)

            mydflistinclude <- mydflist %>% 
            lapply(function(mydf,deltapsicutoff.=deltapsicutoff,containnoveljunction.=containnoveljunction) {
                mydf %>% with((((deltapsi %>% abs)>=deltapsicutoff.) %>% any)&(clusterContainsNovelJunction %>% any)==containnoveljunction.)
            }) %>% 
            unlist
            if ((mydflistinclude %>% sum)==0) return(NULL)

            outdf <- mydflist[mydflistinclude] %>% 
            (function(dflist) do.call(rbind,dflist))
            
            outdf$StartCoord <- outdf$StartCoord %>% as.numeric
            outdf$EndCoord <- outdf$EndCoord %>% as.numeric
            outdf$strand <- gsub('.*_','',outdf$Cluster)
            outdf$refcoord <- paste0(
                outdf$Chr,':',
                ifelse(
                    outdf$strand=='+',
                    paste0('.*-',outdf$StartCoord),
                    paste0(outdf$EndCoord,'-.*')
                )
            )

            outdf <- outdf %>% 
            merge(
                y=EVENT_INFO %>% subset(COMPLEX=='Alt5') %>% 
                mutate(coordo=ifelse(gsub('.*:','',REF_CO)=='-',gsub('-[0-9]*','-.*',COORD_o),gsub(':[0-9]*-',':.*-',COORD_o))),
                #by.x='genes',by.y='GENE',all.x=TRUE
                by.x='refcoord',by.y='coordo',all.x=TRUE
            )

            #message(outdf$EVENT %>% is.na %>% sum)
            outdf
        })
        if (novela5ss %>% is.null) {
            message(paste0('0 alt. 5\' splice sites'))
            return(NULL)
        }

        outdf <- novela5ss %>% 
        data.frame %>% 
        distinct(Chr,JuncStart=StartCoord,JuncEnd=EndCoord,Junction,renamedCluster,genes,splice_site,anchor,known_junction=(known_junction=='1'),ClusterASType,deltapsi,padj,clusterContainsNovelJunction,EVENT,COORD_o,LE_o=gsub(' *','',LE_o) %>% as.numeric,FULL_CO,COMPLEX,REF_CO,LE_n=gsub(' *','',LE_n) %>% as.numeric,CO_C1,CO_A,CO_C2,Seq_C1,Seq_A,Seq_C2) %>% 
        (function(mydf) mydf[,c('Chr','JuncStart','JuncEnd','Junction','renamedCluster','genes','splice_site','anchor','known_junction','ClusterASType','deltapsi','padj','clusterContainsNovelJunction','EVENT','COORD_o','LE_o','FULL_CO','COMPLEX','REF_CO','LE_n','CO_C1','CO_A','CO_C2','Seq_C1','Seq_A','Seq_C2')]) %>% 
        subset((padj %>% as.numeric)<=padjcutoff) %>% 
        merge(EVENT_CONSERVATION,by.x='EVENT',by.y='EventID',all.x=TRUE)

        message(paste0(outdf %>% with(renamedCluster) %>% unique %>% length,' alt. 5\' splice sites'))
        return(outdf)
}

# match alt. 3' splice site
.match_alt3ss <- function(
    intronjunctionsdf=NULL,ASType='A3SS',
    EVENT_INFO=EVENT_INFO.hg38.tab,EVENT_CONSERVATION=EVENT_CONSERVATION.tab.hg38mm10,
    deltapsicutoff=.01,padjcutoff=.1,containnoveljunction=FALSE
    ) {
        novela3ss <- intronjunctionsdf %>% 
        subset(ClusterASType==ASType) %>% 
        (function(mydf) {
            mydflist <- mydf %>% 
            split(f=mydf$renamedCluster)

            mydflistinclude <- mydflist %>% 
            lapply(function(mydf,deltapsicutoff.=deltapsicutoff,containnoveljunction.=containnoveljunction) {
                mydf %>% with((((deltapsi %>% abs)>=deltapsicutoff.) %>% any)&(clusterContainsNovelJunction %>% any)==containnoveljunction.)
            }) %>% 
            unlist
            if ((mydflistinclude %>% sum)==0) return(NULL)

            outdf <- mydflist[mydflistinclude] %>% 
            (function(dflist) do.call(rbind,dflist))

            outdf$StartCoord <- outdf$StartCoord %>% as.character
            outdf$EndCoord <- outdf$EndCoord %>% as.character
            outdf$strand <- gsub('.*_','',outdf$Cluster)
            outdf$refcoord <- paste0(
                outdf$Chr,':',
                ifelse(
                    outdf$strand=='+',
                    paste0(outdf$EndCoord,'-.*'),
                    paste0('.*-',outdf$StartCoord)
                    #paste0('.*',outdf$EndCoord,'.*-.*','','.*'),
                    #paste0('.*','','.*-.*',outdf$StartCoord,'.*')
                )
            )

            outdf <- outdf %>% 
            merge(
                y=EVENT_INFO %>% subset(COMPLEX=='Alt3') %>% 
                mutate(coordo=ifelse(gsub('.*:','',REF_CO)=='+',gsub('-[0-9]*','-.*',COORD_o),gsub(':[0-9]*-',':.*-',COORD_o))),
                #by.x='genes',by.y='GENE',all.x=TRUE
                by.x='refcoord',by.y='coordo',all.x=TRUE
            )
            
            #message(outdf$EVENT %>% is.na %>% sum)
            outdf
        })
        if (novela3ss %>% is.null) {
            message(paste0('0 alt. 3\' splice sites'))
            return(NULL)
        }

        outdf <- novela3ss %>% data.frame %>% 
        distinct(Chr,JuncStart=StartCoord %>% as.numeric,JuncEnd=EndCoord %>% as.numeric,Junction,renamedCluster,genes,splice_site,anchor,known_junction=known_junction=='1',ClusterASType,deltapsi,padj,clusterContainsNovelJunction,EVENT,COORD_o,LE_o=gsub(' *','',LE_o) %>% as.numeric,FULL_CO,COMPLEX,REF_CO,LE_n=gsub(' *','',LE_n) %>% as.numeric,CO_C1,CO_A,CO_C2,Seq_C1,Seq_A,Seq_C2) %>% 
        (function(mydf) mydf[,c('Chr','JuncStart','JuncEnd','Junction','renamedCluster','genes','splice_site','anchor','known_junction','ClusterASType','deltapsi','padj','clusterContainsNovelJunction','EVENT','COORD_o','LE_o','FULL_CO','COMPLEX','REF_CO','LE_n','CO_C1','CO_A','CO_C2','Seq_C1','Seq_A','Seq_C2')]) %>% 
        subset((padj %>% as.numeric)<=padjcutoff) %>% 
        merge(EVENT_CONSERVATION,by.x='EVENT',by.y='EventID',all.x=TRUE)

        message(paste0(outdf %>% with(renamedCluster) %>% unique %>% length,' alt. 3\' splice sites'))
        return(outdf)
}

# match retained intron
.match_retainedintron <- function(
    intronjunctionsdf=NULL,ASType='PossibleRI',
    EVENT_INFO=EVENT_INFO.hg38.tab,EVENT_CONSERVATION=EVENT_CONSERVATION.tab.hg38mm10,
    deltapsicutoff=.01,padjcutoff=.1,containnoveljunction=FALSE
    ) {
        novelretainedintron <- intronjunctionsdf %>% 
        subset(ClusterASType==ASType) %>% 
        (function(mydf) {
            mydflist <- mydf %>% 
            split(f=mydf$renamedCluster)

            mydflistinclude <- mydflist %>% 
            lapply(function(mydf,deltapsicutoff.=deltapsicutoff,containnoveljunction.=containnoveljunction) {
                mydf %>% with((((deltapsi %>% abs)>=deltapsicutoff.) %>% any)&(clusterContainsNovelJunction %>% any)==containnoveljunction.)
            }) %>% 
            unlist
            if ((mydflistinclude %>% sum)==0) return(NULL)

            outdf <- mydflist[mydflistinclude] %>% 
            lapply(function(x) {
                coordfull <- paste0(x$Chr[1],':',x %>% subset(Start==1&End==2) %>% with(StartCoord+1),'-',x %>% subset(Start==1&End==2) %>% with(EndCoord-1))
                coordstartend <- ifelse(
                    test=gsub('.*_','',x$Cluster[1])=='+',
                    yes=paste0(x %>% subset(Start==1&End==2) %>% with(StartCoord),'=',x %>% subset(Start==1&End==2) %>% with(EndCoord)),
                    no=paste0(x %>% subset(Start==1&End==2) %>% with(EndCoord),'=',x %>% subset(Start==1&End==2) %>% with(StartCoord))
                )
                x$coordfull <- gsub(' *','',coordfull)
                x$coordstartend <- gsub(' *','',coordstartend)
                x
            }) %>% 
            (function(dflist) do.call(rbind,dflist))

            outdf <- outdf %>% 
            merge(
                y=EVENT_INFO %>% subset(COMPLEX=='IR'),
                #by.x='genes',by.y='GENE',all.x=TRUE
                by.x='coordfull',by.y='COORD_o',all.x=TRUE
            )
            
            outdf
        })
        if (novelretainedintron %>% is.null) {
            message(paste0('0 retained introns'))
            return(NULL)
        }

        outdf <- novelretainedintron %>% 
        distinct(Chr,JuncStart=StartCoord,JuncEnd=EndCoord,Junction,renamedCluster,genes,splice_site,anchor,known_junction=known_junction=='1',ClusterASType,deltapsi,padj,clusterContainsNovelJunction,EVENT,COORD_o=coordfull,LE_o=gsub(' *','',LE_o) %>% as.numeric,FULL_CO,COMPLEX,REF_CO,LE_n=gsub(' *','',LE_n) %>% as.numeric,CO_C1,CO_A,CO_C2,Seq_C1,Seq_A,Seq_C2) %>% 
        (function(mydf) mydf[,c('Chr','JuncStart','JuncEnd','Junction','renamedCluster','genes','splice_site','anchor','known_junction','ClusterASType','deltapsi','padj','clusterContainsNovelJunction','EVENT','COORD_o','LE_o','FULL_CO','COMPLEX','REF_CO','LE_n','CO_C1','CO_A','CO_C2','Seq_C1','Seq_A','Seq_C2')]) %>% 
        subset((padj %>% as.numeric)<=padjcutoff) %>% 
        merge(EVENT_CONSERVATION,by.x='EVENT',by.y='EventID',all.x=TRUE)  

        message(paste0(outdf %>% with(renamedCluster) %>% unique %>% length,' retained introns'))
        return(outdf)
}