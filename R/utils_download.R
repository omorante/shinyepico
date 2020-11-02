# DOWNLOAD HANDLER FUNCTIONS
fwrite_bed = function(bed_file, file_name, DMR = FALSE, ...) {
  if (!DMR) {
    bed_file = data.table::data.table(
      chr = bed_file$chr,
      start = bed_file$pos - 1,
      end = bed_file$pos,
      name = bed_file$cpg
    )
  }
  
  #We remove NA values in final beds, and remove possible trailing spaces
  bed_file = stats::na.omit(bed_file)
  bed_file$start = trimws(format(as.numeric(bed_file$start), scientific = FALSE))
  bed_file$end = trimws(format(as.numeric(bed_file$end), scientific = FALSE))
  
  data.table::fwrite(
    bed_file,
    file = file_name,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE,
    ...
  )
}


create_dmrs_bed_background = function(mcsea_result,
                                      collapse = FALSE,
                                      regionsTypes,
                                      annotation,
                                      directory){
  
  annotation$cpg = row.names(annotation)
  annotation = data.table::as.data.table(annotation)
  data.table::setkeyv(annotation, "cpg")
  associationTypes = paste0(regionsTypes, "_association")
  
  
  for (association in associationTypes) {
    temp = data.table::as.data.table(t(vapply(names(mcsea_result[[1]][[association]]), function(name) {
      target = annotation[list(mcsea_result[[1]][[association]][[name]]), nomatch = NULL, mult = "all"]
      chr = unique(as.character(target$chr))
      ini = min(as.numeric(target$pos))
      fin = max(as.numeric(target$pos))
      name = paste("background", association, name, sep = "_")
      
      res = c(
        chr = chr,
        start = ini,
        end = fin,
        name = name
      )
      
      if (length(res) != 4) {
        message(
          paste(
            name,
            "DMR is not correctly annotated. NAs will be introduced in final bed"
          )
        )
        res = c(
          chr = "NA",
          start = "NA",
          end = "NA",
          name = name
        )
      }
      
      res
      
    }, character(4))))
    
    if(collapse){
      fwrite_bed(
        temp,
        file_name = paste0(directory, "/", "background.bed"),
        DMR = TRUE, 
        append = TRUE
      )
      
    }
    else{
      fwrite_bed(
        temp,
        file_name = paste0(directory, "/", "background", "_", association, ".bed"),
        DMR = TRUE
      )
      
    }
  }
  
  
}

create_filtered_beds_dmrs = function(mcsea_filtered,
                                     regionsTypes,
                                     annotation,
                                     directory) {
  annotation$cpg = row.names(annotation)
  annotation = data.table::as.data.table(annotation)
  data.table::setkeyv(annotation, "cpg")
  
  associationTypes = paste0(regionsTypes, "_association")
  
  for (contrast in names(mcsea_filtered)) {
    for (association in associationTypes) {
      temp = data.table::as.data.table(t(vapply(names(mcsea_filtered[[contrast]][[association]]), function(name) {
        target = annotation[list(mcsea_filtered[[contrast]][[association]][[name]]), nomatch = NULL, mult = "all"]
        chr = unique(as.character(target$chr))
        ini = min(as.numeric(target$pos))
        fin = max(as.numeric(target$pos))
        name = paste(contrast, association, name, sep = "_")
        
        res = c(
          chr = chr,
          start = ini,
          end = fin,
          name = name
        )
        
        if (length(res) != 4) {
          message(
            paste(
              name,
              "DMR is not correctly annotated. NAs will be introduced in final bed"
            )
          )
          res = c(
            chr = "NA",
            start = "NA",
            end = "NA",
            name = name
          )
        }
        
        res
        
      }, character(4))))
      
      fwrite_bed(
        temp,
        file_name = paste0(directory, "/", contrast, "_", association, ".bed"),
        DMR = TRUE
      )
      
    }
  }
}

create_filtered_bed_dmrs_clusters = function(dendro_data,
                                             mcsea_filtered,
                                             annotation,
                                             directory) {
  annotation$cpg = row.names(annotation)
  annotation = data.table::as.data.table(annotation)
  data.table::setkeyv(annotation, "cpg")
  
  lapply(unique(dendro_data), function(cluster) {
    temp = data.table::as.data.table(t(vapply(names(dendro_data[dendro_data == cluster]), function(name) {
      contrast = limma::strsplit2(name, "\\|")[1]
      association = paste0(limma::strsplit2(name, "\\|")[2], "_association")
      gene = limma::strsplit2(name, "\\|")[3]
      
      target = annotation[list(mcsea_filtered[[contrast]][[association]][[gene]]), nomatch = NULL, mult = "all"]
      
      chr = unique(as.character(target$chr))
      ini = min(as.numeric(target$pos))
      fin = max(as.numeric(target$pos))
      name = paste(contrast, association, gene, sep = "_")
      
      res = c(
        chr = chr,
        start = ini,
        end = fin,
        name = name
      )
      
      if (length(res) != 4) {
        message(paste(
          name,
          "DMR is not correctly annotated. NAs will be introduced in final bed"
        ))
        res = c(
          chr = "NA",
          start = "NA",
          end = "NA",
          name = name
        )
      }
      
      res
      
    }, character(4))))
    
    fwrite_bed(
      temp,
      file_name = paste0(directory, "/", "cluster", "_", which(unique(dendro_data) %in% cluster), ".bed"),
      DMR = TRUE
    )
    
  }
  
  
  
  )
  
  
  
}

create_filtered_beds = function(filtered_data, annotation, directory) {
  
  annotation$cpg = row.names(annotation)
  annotation = data.table::as.data.table(annotation)
  
  #saving hypo and hyper results individually
  lapply(names(filtered_data), function(name) {
    temp = data.table::merge.data.table(filtered_data[[name]][filtered_data[[name]]$dif_current < 0,],
                                        annotation,
                                        by = "cpg",
                                        all.x = TRUE)
    
    fwrite_bed(temp,
               file_name = paste0(directory, "/", name, "_hypermethylated.bed"))
  })
  
  lapply(names(filtered_data), function(name) {
    temp = data.table::merge.data.table(filtered_data[[name]][filtered_data[[name]]$dif_current > 0, ],
                                        annotation,
                                        by = "cpg",
                                        all.x = TRUE)
    
    fwrite_bed(temp,
               file_name = paste0(directory, "/", name, "_hypomethylated.bed"))
  })
  
  
  fwrite_bed(annotation, file_name = paste0(directory, "/", "annotation.bed"))
  
}

create_filtered_bed_clusters = function(dendro_data, annotation, directory) {
  annotation$cpg = row.names(annotation)
  annotation = data.table::setDT(annotation)
  
  #saving results by cluster
  lapply(unique(dendro_data), function(cluster) {
    temp = annotation[annotation$cpg %in% names(dendro_data)[dendro_data == cluster], ]
    
    fwrite_bed(temp,
               file_name = paste0(
                 directory,
                 "/",
                 "Cluster_",
                 which(unique(dendro_data) %in% cluster),
                 ".bed"
               ))
  })
  
  fwrite_bed(annotation, file_name = paste0(directory, "/", "annotation.bed"))
}
