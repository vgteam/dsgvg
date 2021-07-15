.. _odgi tips:

#########
odgi tips
#########

Identifying break point positions relative to given query (reference) path(s) of all the tips in the graph or of tips of given path(s). Prints BED records to stdout.

SYNOPSIS
========

**odgi tips** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi tips command identifies break point positions relative to given query (reference) path(s) of all the tips in
the graph or of tips of given path(s). Prints BED records to stdout. Each record consist of:

- **chrom**: The query path name.
- **start**: The 0-based start position of the query we hit in the node.
- **end**: The 1-based end position of the query we hit in the node.
- **median_range**: The 0-based median of the whole query path range of the node we hit. It is possible that a node contains several steps, so we want to mirror that here.
- **path**: The name of the path we walked.
- **path_pos**: The 0-based position of the path we walked when we hit the node of the query path.
- **path_pos**: The 0-based position of the path we walked when we hit the node of the query path.
- **walk_from_front**: If `1` we walked from the front of the target path. Else it is `0`.

OPTIONS
=======

MANDATORY OPTIONS
-----------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!


Tips Options
-------------

| **-q, --query-path**\ =\ *NAME*
| Use this query path.

| **-r, --target-path**\ =\ *NAME*
| Use this target (reference) path.

| **-Q, --query-paths**\ =\ *FILE*
| Use query paths listed in *FILE*. Each line must contain one path.

| **-R, --target-paths**\ =\ *FILE*
| Use target (reference) paths list (one per line) in *FILE*.

| **-v, --not-visited-tsv**\ =\ *FILE*
| Write target path(s) that do not visit the query path(s) to this *FILE*.

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
| Print a help message for **odgi tips**.

..
	EXIT STATUS
	===========

	| **0**
	| Success.

	| **1**
	| Failure (syntax or usage error; parameter error; file processing
		failure; unexpected error).
..
	BUGS
	====

	Refer to the **odgi** issue tracker at
	https://github.com/pangenome/odgi/issues.
