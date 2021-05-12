# pan-gene stats 

pan_stats <- read.csv(file = "pan_matrix_stats.csv", header= TRUE)
pan_stats


# get summary of number of pan gene per pan gene class 
aggregate(pan_stats$n, by=list(class=pan_stats$class), FUN=sum)


#syntenic section

# core,near core 
core_near_core = pan_stats %>% filter(class=="Core Gene" | class=="Near-Core Gene")
aggregate(core_near_core$n, by=list(Subgenome=core_near_core$Subgenome), FUN=sum)
# Subgenome     x
# 1           M1 11351
# 2           M2  7168
# M1 + M2 18519
# 3 Non-Syntenic 13533
# 18519/(18519+13533) = 0.5777799

dispensable_private = pan_stats %>% filter(class=="Dispensable Gene" | class=="Private Gene")
aggregate(dispensable_private$n, by=list(Subgenome=dispensable_private$Subgenome), FUN=sum)
#     Subgenome     x
# 1           M1   836
# 2           M2   465
# 3 Non-Syntenic 69680
# (836+465)/(836+465+69680) = 0.01832885

aggregate(pan_stats$n, by=list(class=pan_stats$class), FUN=sum)

# total M1M2 in core and near core 18519
# total M1M2 in dispensable and private 1301
# total M1M2 19820

# split by core and near core 
core = pan_stats %>% filter(class=="Core Gene")
aggregate(core$n, by=list(Subgenome=core$Subgenome), FUN=sum)
#      Subgenome     x
# 1           M1 10435
# 2           M2  6857
# 3 Non-Syntenic 10618
# (10435+6857)/19820

near_core = pan_stats %>% filter(class=="Near-Core Gene")
aggregate(near_core$n, by=list(Subgenome=near_core$Subgenome), FUN=sum)

#     Subgenome    x
# 1           M1  916
# 2           M2  311
# 3 Non-Syntenic 2915
# (916+311)/19820
