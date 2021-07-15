.. _odgi paths:

#########
odgi paths
#########

Interrogate the embedded paths of a graph. Does not print anything to stdout by default!

SYNOPSIS
========

**odgi paths** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi paths command allows the investigation of paths of a given
variation graph. It can calculate overlap statistics of groupings of
paths.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Path Investigation Options
---------------------

| **-O, --overlaps**\ =\ *FILE*
| Read in the path grouping *FILE* to generate the overlap statistics
  from. The file must be tab-delimited. The first column lists a
  grouping and the second the path itself. Each line has one path entry.
  For each group the pairwise overlap statistics for each pairing will
  be calculated and printed to stdout.

| **-L, --list-paths**
| Print the paths in the graph to stdout. Each path is printed in its
  own line.

| **-l, --list-paths-start-end**
| If **-L,--list-paths** was specified, this additionally prints the start and end positions of each path in additional, tab-delimited coloumns.

| **-H, --haplotypes**
| Print to stdout the paths in an approximate binary haplotype matrix
  based on the graph’s sort order. The output is tab-delimited:
  **path.name**, *path.length*, *node.count*, *node.1*,
  *node.2*, *node.n*. Each path entry is printed in its own line.

| **-D, --delim**\ =\ *CHAR*
| The part of each path name before this delimiter is a group
  identifier. This parameter should only be set in combination with
  **-H, --haplotypes**. Prints an additional, first column
  **group.name** to stdout.

| **-d, --distance**
| Provides a sparse distance matrix for paths. If **-D, --delim** is
  set, it will be path groups distances. Each line prints in a tab-delimited format to stdout:
  *path.a*, *path.b*, *path.a.length*, *path.b.length*, *intersection*, *jaccard*, *euclidean*.

| **-f, --fasta**
| Print paths in FASTA format to stdout. One line for the FASTA header, another line for the whole sequence.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi paths**.

..
	EXIT STATUS
	===========
	
	| **0**
	| Success.
	
	| **1**
	| Failure (syntax or usage error; parameter error; file processing
	  failure; unexpected error).
	
	BUGS
	====
	
	Refer to the **odgi** issue tracker at
	https://github.com/pangenome/odgi/issues.
