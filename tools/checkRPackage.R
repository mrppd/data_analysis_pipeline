args = commandArgs(trailingOnly=TRUE)

if(length(args)==1){
  pkName = args[1]
  res = pkName %in% rownames(installed.packages())
  if(res==TRUE){
    print(1)
  } else {
    print(0)
  }
} else {
  print(0)
}



