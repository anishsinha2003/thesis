# Make cells comparable. When sequencing happens, the machine does NOT directly measure:
# “this gene has value 10”. Instead, it reads tiny fragments of RNA millions of times
# and each fragment detected is called a "read". More read means more sequencing coverage
# So before normalisation the value at the expression value of each gene in the cells (in the matrix)
# is how many reads mapped to each gene (basically read counts).
# however, different cells often get different sequencing depth (read counts)
# so if cell A has total 50,000 reads and cell B has 5000 reads This DOES NOT necessarily mean: Cell A is biologically more active
# It might simply mean: the sequencer captured more RNA fragments from Cell A
# Without normalisation: Cell A appears more “expressed” even if biologically: Cell A and Cell B are actually similar
# What normalisation does: It rescales expression values so: all cells are on similar scale


