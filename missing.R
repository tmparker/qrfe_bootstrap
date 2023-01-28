# Check the parts directory for missing runs, since sometimes they get dropped
# on the big computer.

nstream <- 200
des <- 1:10
part <- 1:nstream
all_parts <- 1:2000

# Find which jobs have been dropped.
parts <- paste0("part", part, ".rda")
big_list <- c(sapply(des, function(x) paste0("design", x, parts)))
# done <- match(list.files("./parts"), big_list)
done <- match(readLines("miss.txt"), big_list) # from ls parts/ > miss.txt
redo <- all_parts[-done] - 1
# the redo vector tells you which jobs to send back to the system.
# numbering is different from R numbering

dmiss <- ceiling(redo / 200) # which designs are missing


# Print the job numbers for those that remain to be done in a way that can be
# easily pasted into the run_simulations.sh script.
jmp <- which(diff(redo) != 1)
if (jmp[length(jmp)] != redo[length(redo)]) {jmp <- c(jmp, length(redo))}
contig <- vector(mode = "list", length = length(jmp))
contig[[1]] <- seq(redo[1], redo[1] + jmp[1] - 1)
for (i in 2:length(jmp)) {
  contig[[i]] <- seq(redo[jmp[i - 1] + 1], redo[jmp[i]])
}
fun <- function(s) {
  if (length(s) == 1) { paste0(s, ",", collapse = "")
  } else { paste0(s[1], "-", s[length(s)], ",", collapse = "") }
}
all_one <- paste(sapply(contig, fun), sep = "", collapse = "")
no_comma <- substr(all_one, 1, nchar(all_one) - 1)
cat(no_comma)

