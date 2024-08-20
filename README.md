# PangraphViz

This is a visualization for data produced by <a href="https://neherlab.github.io/pangraph/tutorials/tutorial_1/">pangraph</a>, which you can read more about in these <a href="https://pubmed.ncbi.nlm.nih.gov/37278719/">two</a> <a href="https://www.biorxiv.org/content/10.1101/2024.07.08.602537v1.full.pdf">papers</a> from Nicholas Noll, Marco Molari, Liam Shaw, and Richard Neher. A pangraph is a way to look at a population of closely related genomes, and this browser is an attempt to make that <i>looking</i> easier. In brief, a pangraph represents each genome as some path through a series of blocks of DNA, and each block represents a multiple sequence alignment of very closely related sequences. Some of these blocks appear only once and appear in every genome - these are called "core" blocks. Other blocks are either not present in every genome or present in every genome but repeated in multiple locations.

This visualization was made by Milo Johnson (milo.s.johnson.13@gmail.com) as part of a group project at the Kavli Institute for Theoretical Physics with help from Marco Molari, Richard Neher, Indra Gonzalez Ojeda, Caelan Brooks, Sofya Garushyants, Isaac Gifford, Athina Diakogianni, and James Ferrare.

There are very likely to be bugs here, so feel free to leave issues or email about them or about ideas for how this could be made more useful. The current to-do / to-add list is:

TODO 
* make visual gaps in core genome for junctions using local pangenome size
* make an alignment viewer for when you click on a block
* work on annotation (quality, merging, pipeline, etc.)
* add original e coli dataset from paper
* general code clean up