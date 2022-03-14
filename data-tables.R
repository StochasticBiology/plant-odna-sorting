# read heteroplasmy observations from Excel files and parse for use in analysis

library(readxl)

# read data from Excel sheets and covert to dataframes
# between-generation data
gen.tib = read_excel("MSH1_heteroplasmy_progeny_22.xlsx")
gen.df = data.frame(gen.tib)
# within-child data
tissue.tib = read_excel("Within_plant_ALL_data_22.xlsx")
tissue.df = data.frame(tissue.tib)

# loop through organelle-background sets
organelles = c("mito", "plastid")
backgrounds = c("MSH1", "WILD")

# function to deal with awkward measurements (presumably due to experimental noise)
h.normalise = function(h) {
  if(h > 99.5) { return(1) }
  if(h < 0.5) { return(0) }
  return(h/100)
}

# extract useful information for between-generation data
gen.set = subset(gen.df, select=c(organelle, unique.familiy.ID, individual_.altSNV, background))
colnames(gen.set) <- c("organelle", "family", "h", "background")

# Table 2: build two (MT+PT) lists of heteroplasmy samples for MSH1 background, between generation data
t2.df = gen.set[gen.set$background == "MSH1" & !is.na(gen.set$background),]
t2.mt = t2.pt = list()
for(family in unique(t2.df$family)) {
  sub = t2.df[t2.df$family == family,]
  if(sub$organelle[1] == "mito") {
    t2.mt[[length(t2.mt)+1]] = unlist(lapply(sub$h, h.normalise))
  } else {
    t2.pt[[length(t2.pt)+1]] = unlist(lapply(sub$h, h.normalise))
  }
}

# Table 4: build one (MT) list of heteroplasmy samples for WT background, between generation data
t4.df = gen.set[gen.set$background == "WILD" & !is.na(gen.set$background),]
t4.mt = list()
for(family in unique(t4.df$family)) {
  sub = t4.df[t4.df$family == family,]
  t4.mt[[length(t4.mt)+1]] = unlist(lapply(sub$h, h.normalise))
}

# extract useful information for within-individual data
tissue.set = subset(tissue.df, select=c(organelle, uniqueptID, X.alt_corrected, background))
colnames(tissue.set) <- c("organelle", "individual", "h", "background")

# Table 3: build two (MT+PT) list of heteroplasmy samples for MSH1 background, within-individual data
t3.df = tissue.set[!is.na(tissue.set$background),]
t3.mt = t3.pt = list()
for(individual in unique(t3.df$individual)) {
  sub = t3.df[t3.df$individual == individual,]
  if(sub$organelle[1] == "mito") {
    t3.mt[[length(t3.mt)+1]] = unlist(lapply(sub$h, h.normalise))
  } else {
    t3.pt[[length(t3.pt)+1]] = unlist(lapply(sub$h, h.normalise))
  }
}

