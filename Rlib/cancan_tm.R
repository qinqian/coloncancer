################################
## text mining on cancer
################################

library(tm)
txt <- system.file("texts", "txt", package = "tm")
(ovid <- Corpus(DirSource(txt, encoding = "UTF-8"),
                readerControl = list(language="lat")))

docs <- c("This is a text.", "This is another one")

Corpus(VectorSource(docs))
reut <- system.file("texts", "crude", package = "tm")

reuter <- Corpus(DirSource(reut),
                 readerControl = list(reader = readReut21578XML))

writeCorpus(ovid, path="test")

inspect(ovid[1:2])

identical(ovid[[2]], ovid[["ovid_2.txt"]])
reuters <- tm_map(reuter, as.PlainTextDocument)

reuters <- tm_map(reuters, tolower)
reuters <- tm_map(reuters, stripWhitespace)


reuters <- tm_map(reuters, removeWords, stopwords("english"))

library(rJava)
library(Snowball)
tm_map(reuters, stemDocument)

query <- "id == '237' & heading == 'INDONESIA SEEN AT CROSSROADS OVER ECONOMIC CHANGE'"

tm_filter(reuters, FUN = sFilter, query)

library(meta) ## need install
library(metafor)

data(crude)
DublinCore(crude[[1]], "Creator") <- "Ano Nymous"
meta(crude[[1]])

meta(crude, tag = "test", type = "corpus") <- "test meta"
meta(crude, type = "corpus")

meta(crude, "foo") <- letters[1:20]
meta(crude)

dtm <- DocumentTermMatrix(reuters)
inspect(dtm[1:5,100:105])
