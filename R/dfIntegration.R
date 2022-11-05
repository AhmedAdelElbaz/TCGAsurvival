#' Integrate Data Frames
integrate_files <- function(to_dataframe  ,
                            from_dataframe ,
                            key_to ,
                            key_from ,
                            columns.to.check=key_from,
                            empty.values = 'NULL')
{
  library(hash)
  to_dataframe <- as.data.frame(to_dataframe)
  from_dataframe <- as.data.frame(from_dataframe)
  for (each.column in columns.to.check)
  {
    to_dataframe[each.column] = empty.values
  }

  #Create a dictionary of the dataframe we wanna copy values from
  dict_from = hash()
  for (i in 1:nrow(from_dataframe)) #i = each row
  {
    #list_of_values = list()
    list_of_values = hash()
    for (each.column in columns.to.check)
    {
      list_of_values[[each.column]] = from_dataframe[i,each.column]
    }

    dict_from[[from_dataframe[i,key_from]]] = list_of_values
  }

  for (i in 1:nrow(to_dataframe))
  {
    key = to_dataframe[i,key_to]
    for (each.column in columns.to.check)
    {
      to_dataframe[i,each.column] = dict_from[[key]][[each.column]]
    }
  }
  return(to_dataframe)
}
