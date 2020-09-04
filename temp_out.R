test = rbind(table(df.lasR_meta$TRUNCATED[df.lasR_meta$GROUP == 'CF']), table(df.lasR_meta$TRUNCATED[df.lasR_meta$GROUP == 'ENV']), table(df.lasR_meta$TRUNCATED[df.lasR_meta$GROUP == 'WND']))
rownames(test) = c('CF', 'ENV', 'WND')
test = as.data.frame(test)

test$nIsolates = test$FULL + test$TRUNC
test$pct_FULL = test$FULL/test$nIsolates
test$pct_TRUNC = test$TRUNC/test$nIsolates