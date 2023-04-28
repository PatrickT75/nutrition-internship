# The function get_description takes three arguments: data file, reference file, and name of new file.
# First row of data file will be skipped.
get_descriptions = function(datafile, referencefile, newfile) {
  data = read.delim(datafile, header = TRUE, skip=1)
  
  ref = read.delim(referencefile, header = FALSE)
  colnames(ref) = c('Pathway', 'Description')
  refpwy = ref["Pathway"]
  refdesc = ref["Description"]
  
  pwy_names = data[,1]
  pwy_descriptions = rep('NA', length(pwy_names))
  
  a = 1
  for (name in pwy_names) {
    ind = which(refpwy == name)
    desc = refdesc[ind,]
    pwy_descriptions[a] = desc
    a = a+1
  }
  
  data_pwy = data[1]
  data_rem = data[2:ncol(data)]
  data_desc = cbind(data_pwy, Description = pwy_descriptions)
  data_desc = cbind(data_desc, data_rem)
  
  write.table(data_desc, newfile, sep = '\t', row.names = FALSE)
}