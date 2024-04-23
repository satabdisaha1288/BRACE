#Split into train test for data replicates
split_train_test<- function(split_ratio, data){
  # Sample dataset (replace this with your actual dataset)
  set.seed(123)
  num_rows <- nrow(data)
  num_train <- round(num_rows * split_ratio)
  num_test <- num_rows - num_train
  # Create random indices for the training and testing sets
  train_indices <- sample(1:num_rows, num_train, replace = FALSE)
  test_indices <- setdiff(1:num_rows, train_indices)
  # Split the data into training and testing sets
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  return(list(train_data = train_data, test_data = test_data))
}

#Split all replicates in train test
split_replicates_test_train<- function(split_ratio, XlogRelAbun, XRelAbun, y){
  y_split = lapply(y, function(x) split_train_test(split_ratio = split_ratio , data = x))
  XlogRelAbun_split = lapply(XlogRelAbun, function(x) 
    split_train_test(split_ratio = split_ratio , data = x))
  XRelAbun_split = lapply(XRelAbun, function(x) 
    split_train_test(split_ratio = split_ratio , data = x))
  # Xabs_split = lapply(Xabs, function(x) 
  #   split_train_test(split_ratio = split_ratio , data = x))
  return(list(y_split = y_split, XlogRelAbun_split = XlogRelAbun_split, 
         XRelAbun_split = XRelAbun_split))
}
