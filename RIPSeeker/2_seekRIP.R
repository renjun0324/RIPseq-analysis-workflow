
chr = paste0("chr", c(1:22,"M","X","Y"))
binsize = 200
treat = "HA-IP"
control = "IgG-IP"
treat_path = paste0(root_path, "/7_RIPSeeker/", treat, "/binsize", binsize)
control_path = paste0(root_path, "/7_RIPSeeker/", control, "/binsize", binsize, "/")

save_path = paste0(root_path, "/7_RIPSeeker/result_", treat, "_", control, "/binsize", binsize, "/chr")
system(paste0("mkdir -p ", save_path))

#-----------------------------------------------------------------
#                                                                
#                       parameter setting
#                                                                 
#-----------------------------------------------------------------

## seekRIP
padjMethod="BH"
logOddCutoff=-Inf
pvalCutoff=1
pvalAdjCutoff=1
eFDRCutoff=1

#-----------------------------------------------------------------
#                                                                
#                            seekrip
#                                                                 
#-----------------------------------------------------------------

for(i in 1:25){
  cat(chr[i], " ")
  
  treat_rda = paste0(treat, "_", chr[i])
  control_rda = paste0(control, "_", chr[i])
  cat(treat_rda, "|", control_rda, "\n")
  
  load(paste0(treat_path, "/", chr[i],"/",treat_rda,".rda"))
  load(paste0(control_path, "/", chr[i],"/",control_rda,".rda"))
  
  ## running
  cat(median(width(get(treat_rda)$nbhGRList)), "\n")
  seekrip = seekRIP(nbhGRRIP = get(treat_rda)$nbhGRList,
                    nbhGRCTL = get(control_rda)$nbhGRList)
  
  ## save
  name = paste0("seekrip_",chr[i])
  assign(name, seekrip)
  save(list = name, file = paste0(save_path, "/", name,".rda"))
}

#-----------------------------------------------------------------
#                                                                
#                            merge
#                                                                 
#-----------------------------------------------------------------

riplist <- lapply(1:25, function(i){
  seekrip = paste0("seekrip_",chr[i])
  cat(seekrip, "\n")
  load(paste0(save_path, "/", seekrip,".rda"))
  as.data.frame(get(seekrip))
})
rip <- do.call(rbind, riplist) %>% as.data.frame
ripGR <- as(rip,"GRanges")
save(ripGR, file = paste0(save_path,"ripGR.rda"))
