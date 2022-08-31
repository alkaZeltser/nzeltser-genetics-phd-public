##### format-cancer-nordic-twin-heritability-data.R ##################

# This script parses heritability data published by the nordic twin cancer study
# and plots heritability estimates by cancer type in a forest plot.

library(BoutrosLab.plotting.general);

TABLE <- '/Users/nzeltser/Desktop/Boutros/Germline-somatic/tabula-nordic_twin_cancer.csv';
OUTPUT <- '/Users/nzeltser/Coding-sandbox/germline-somatic/plots/'

# function that parses original value in norcan heritability - familial environment column
parse.norcan.heritability <- function(string) {
  components <- unlist(strsplit(string, split = ' '));
  
  norcan.df <- data.frame(
    norcan.heritability = components[1],
    norcan.h.CI = components[2],
    norcan.shared.environment = components[3],
    norcan.shrdenv.CI = components[4]
    );
  
  return(norcan.df);
  
}

# function that parses the original confidence interval string into
# numeric components (left and right confidence interval length).
split.confidence.interval <- function(ci, mean) {
  
  components <- as.numeric(unlist(strsplit(ci, split = '-')));
  
  split.interval <- data.frame(
    left.length = as.numeric(mean) - as.numeric(components[1]),
    right.length = as.numeric(components[2]) - as.numeric(mean)
  );
  
  return(split.interval);
  
}

table <- read.table(TABLE, header = FALSE, sep = ',', skip = 5);
colnames(table) <- c('cancer', 'heritability.familial', 'lichtenstein.heritability', 'n.mono.con', 'm.mono.dis', 'n.di.con', 'n.di.dis', 'sample.size');

# separate heritability from shared environment from confidence intervals
parsed.table <- do.call(rbind, lapply(table$heritability.familial, FUN = parse.norcan.heritability))

# remove percentage sign characters
parsed.table <- data.frame(apply(parsed.table, MARGIN = 2, FUN = function(x) { gsub('%', '', x)}))

# remove parenthesis characters
parsed.table <- data.frame(apply(parsed.table, MARGIN = 2, FUN = function(x) { gsub('\\(', '', x)}))
parsed.table <- data.frame(apply(parsed.table, MARGIN = 2, FUN = function(x) { gsub('\\)', '', x)}))

# split confidence interval ranges into length to the right and left of the center
heritability.split.ci <- do.call(rbind, apply(parsed.table, MARGIN = 1, FUN = function(x) {split.confidence.interval(ci = x[2], mean = x[1])}))
colnames(heritability.split.ci) <- paste0('herit.ci.', colnames(heritability.split.ci))
familial.env.split.ci <- do.call(rbind, apply(parsed.table, MARGIN = 1, FUN = function(x) {split.confidence.interval(ci = x[4], mean = x[3])}))
colnames(familial.env.split.ci) <- paste0('shrdenv.ci.', colnames(familial.env.split.ci))

# combine into one df
parsed.table <- data.frame(parsed.table, heritability.split.ci, familial.env.split.ci)

# convert variables to numeric
parsed.table$norcan.heritability <- as.numeric(parsed.table$norcan.heritability);

# recombine parsed data with cancer name category
parsed.table$cancer <- table$cancer;

# recombine parsed data with sample size
parsed.table$sample.size <- table$sample.size;

# separate out overall cancer
overall.cancer <- subset(parsed.table, cancer == 'Overall cancer');
individual.cancer <- subset(parsed.table, cancer != 'Overall cancer');

# order by heritability point estimate
individual.cancer <- individual.cancer[order(individual.cancer$norcan.heritability), ];

# add placeholder number for each cancer to determine plotting order.
# leave space for 'overall cancer' to be plotted last (placeholder == 1)
individual.cancer$placeholder <- c(1:nrow(individual.cancer)) + 1;

overall.cancer$placeholder <- 1;

#bind overall.cancer back into main df
all.cancer <- rbind(overall.cancer, individual.cancer);


# plot forest plot for heritability
create.scatterplot(
  placeholder ~ as.numeric(norcan.heritability),
  data = all.cancer,
  main = 'Mucci et al., NorTwinCan heritability',
  main.cex = 1.5,
#  filename = paste0(OUTPUT, 'cancer-heritability-forest.tiff'),
  x.error.right = all.cancer$herit.ci.right.length,
  x.error.left = all.cancer$herit.ci.left.length,
  error.bar.length = 0,
  xlimits = c(-5, 150),
  ylimits = c(0, max(all.cancer$placeholder) + 2),
  yat = all.cancer$placeholder,
  xat = seq(0, 100, 20),
  yaxis.lab = all.cancer$cancer,
  yaxis.cex = 1,
  ylab.label = 'cancer',
  xlab.label = '% heritability',
  width = 6,
  add.text = TRUE,
  text.x = rep(120, nrow(all.cancer) + 1),
  text.y = c(all.cancer$placeholder, max(all.cancer$placeholder) + 1),
  text.labels = c(all.cancer$sample.size, 'n'),
  text.cex = 0.75,
  text.fontface = c(rep('plain', nrow(all.cancer)), 'bold')
);
plot(plot.lol)
View(plot.lol)


