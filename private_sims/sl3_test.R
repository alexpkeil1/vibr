library(vibr)
library(sl3)

data(metals, package="qgcomp")
set.seed(NULL)
(task1 <- sl3_Task$new(data=metals, covariates=names(metals)[1:23], outcome="y", folds=10))
(task2 <- sl3_Task$new(data=metals, covariates=names(metals)[1:23], outcome="y", folds=task1$folds))


rbind(task1$folds[[1]]$validation_set,
      task2$folds[[1]]$validation_set)

maketassks <- function(){
  currseed = .Random.seed
  set.seed(currseed[1])
  (task1 <- sl3_Task$new(data=metals, covariates=names(metals)[1:23], outcome="y", folds=10))$folds
  (task2 <- sl3_Task$new(data=metals, covariates=names(metals)[1:23], outcome="y", folds=task1$folds))
  .Random.seed = currseed
  list(task1, task2)
}


rm(".Random.seed")
set.seed(NULL)
tsk <- maketassks()

rbind(tsk[[1]]$folds[[1]]$validation_set,
      tsk[[2]]$folds[[1]]$validation_set)
