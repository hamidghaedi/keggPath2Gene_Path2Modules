# keggPath2Gene_Path2Modules
R script to map gene name to each KEGG pathways and Modules. This code can help to have gene name in one column and KEGG pathway ID in another column (```gen2path```). 
Another data set that generates is genes and their  KEGG modules (```gene2module```)
Several intermediate dataset also will be generated during the process that might be of interest, like a dataset on modules, their definition, and pathways involved in that module (```keggModule```).

```R
    library(org.Hs.eg.db)
    library(dplyr)
    library(tidyr)
    library(jsonlite)
    
    db <- org.Hs.egPATH
    # Get the entrez gene identifiers that are mapped to a KEGG pathway ID
    mapped_genes <- mappedkeys(db)
    # Convert to a list
    mapped_genesList <- as.list(db[mapped_genes])
    # converting list to dataframe
    df <- plyr::ldply (mapped_genesList, data.frame)
    mappedDF = data.frame(ENTREZ_ID = as.numeric(df[,1]), KEGG_ID = paste0("map",df[,2]))
    # gen2pathway dataset
     gen2path = aggregate(. ~ ENTREZ_ID, mappedDF, FUN = function(x) 
      toString(x), na.action = NULL)

    # Retriving KEGG module data
    url = "https://www.genome.jp/kegg-bin/download_htext?htext=ko00002&format=json&filedir="
    download.file(url, destfile = "~/keggM.json", method = "curl")
    
    # reading json
    document <- fromJSON(txt=url)
    # parsing json
    df = data.frame(Reduce(rbind, document))
    # pathway modules
    pathMod = df[2,2]
    pathModDF = data.frame(Reduce(rbind, pathMod))
    pathway_modules = data.frame(name =c(), modules = c(), path = c(), p1Path = c(), p2Path = c())
    for(f in 1:dim(pathModDF)[1]){
      for(i in 1:dim(pathModDF)[1]){
        tmp = pathModDF[[2]][[i]]
        for(j in 1:dim(tmp)[1]){
          tmp2 = tmp[[2]][[j]]
          tmp2$module = substr(tmp2$name,1,7)
          tmp2$path = stringr::str_extract(string = tmp2$name, pattern = "(?<=\\[)[^{}]+(?=\\])")
          tmp2$path = sub("PATH:", "", tmp2$path)
          tmp2$p1Path = tmp[j,1]
          tmp2$p2Path = pathModDF[f,1]
          pathway_modules= rbind(pathway_modules, tmp2)
        }
      }
    }
    
    sigMod = df[3,2]
    sigModDF = data.frame(Reduce(rbind, sigMod))
    sig_modules = data.frame(name =c(), modules = c(), path = c(), p1Path = c(), p2Path = c())
    for(f in 1:dim(sigModDF)[1]){
      for(i in 1:dim(sigModDF)[1]){
        tmp = sigModDF[[2]][[i]]
        for(j in 1:dim(tmp)[1]){
          tmp2 = tmp[[2]][[j]]
          tmp2$module = substr(tmp2$name,1,7)
          tmp2$path = stringr::str_extract(string = tmp2$name, pattern = "(?<=\\[)[^{}]+(?=\\])")
          tmp2$path = sub("PATH:", "", tmp2$path)
          tmp2$p1Path = tmp[j,1]
          tmp2$p2Path = sigModDF[f,1]
          sig_modules= rbind(sig_modules, tmp2)
          
        }
      }
    }
    
    pathway_modules$module_type = "pathway"
    sig_modules$module_type = "signature"
    
    keggModule = rbind(pathway_modules, sig_modules) 
    
    # module matrix
    modMat = keggModule[,c(2,3)]
    modMat = data.frame(cbind(keggModule[,2], stringr::str_split_fixed(keggModule$path, " ", 7))) # 7 maximum pathways assigned to a Module
    modMat[modMat == ""] <- NA
    # long dataframe fro module and pathway
    path2Mod = data.frame(module = c(), path = c())
    for(i in 2:ncol(modMat)){
      tmp = modMat[,c(1,i)][!is.na(modMat[i]),]
      names(tmp) = c("module", "path")
      path2Mod = rbind(path2Mod, tmp)
    }
    
    #deduplication
    path2Mod= path2Mod[!duplicated(paste0(path2Mod$module, path2Mod$path)),]
    
    # joing datasets
    kegg = merge(path2Mod, mappedDF, by.x = "path" , by.y ="KEGG_ID" , all.x = TRUE)
    # removing duplicates in kegg dataset
    keggdedup = kegg[!duplicated(paste0(kegg$module, kegg$ENTREZ_ID)),]
    
    #gene2module dataset
    gene2module = aggregate(. ~ ENTREZ_ID, keggdedup[, c(2,3)], FUN = function(x) 
      toString(x), na.action = NULL)
     
     ```
