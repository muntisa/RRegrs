
#assume the data frame below is the training set, with first column being
#net cell association (y), which is the same assumption throughtout RRgres code

my.datf.train

hat(my.datf.train[,2:dim(my.datf.train)[2]])#gives you the hat values!


#########
#########
#IF model is lm, glm, then the following commend dies the same thing
#(again i'm using notation from the original code i've uploaded)
hatvalues(lm.fit$finalModel)

#OTHER things that can be calculated are, e.g.

cooks.distance(lm.fit$finalModel)
#which is the sattistic plotted in standard lm scatterplots

influence(lm.fit$finalModel)#produces several statistics of the kind


#########
#########
#PACKAGE 'ouliers' gives you options for outliers tests which is pretty good
#


########
########
#THRESHOLD values: 3m/n
#where m is the number of parameters, and n number of observations
#

thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1]

