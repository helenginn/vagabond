\page limits List of known limitations

This page contains the list of known limitations (and then a list of known 
issues or bugs) for which a plan to remedy the situation is not yet established.
This is mainly due to being limited in womanpower.

## Known limitations

* Not all PDBs are perfectly imported into Vagabond due to various special 
cases.
* Vagabond is not able to handle nucleic acids.
* Vagabond cannot handle custom geometry for individual ligands.
* Vagabond does not handle real space maps such as those created by cryo-
electron microscopy.
* Vagabond has only been tested properly on Linux, is mildly unstable on Mac OS X
and will probably fall apart on Windows.

## Known bugs

* Vagabond doesn't perfectly apply symmetry operations to certain space groups.
* Post-translational modifications aren't imported and cause a chain to be
dropped if they are listed as HETATM entries (fixable by hand).
* Prolines are not given different geometry if they are cis/trans.

Return to \ref index "the main page".
