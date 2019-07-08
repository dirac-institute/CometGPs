The binary asteroid data are stored in four files:

binastd_PUB.txt contains best estimates of the compiled parameters.

binaste_PUB.txt gives uncertainties of the estimates given in binast_d.txt.
             Either a "+/-" or a "factor of" uncertainty (marked with "*" preceding the error value) is given.
             When an upper limit on "+/-" uncertainty is given, it is marked with the "<" sign.

binastr_PUB.txt gives references and notes for the compiled estimates.

binastm_PUB.txt gives a quick info for each of the estimates on whether it has been estimated directly from observations
             ("measu"), or whether it has been derived from other estimated values ("deriv"; for details on such items,
             see the notes and the paper), or if it has been assumed ("assum"; see comments in the paper).
             In some cases where less reliable data were obtained (often with some unusual assumptions made in their estimations in the original publications),
             the method remark is given in parentheses: "(mea)".  When a derived value was affected by such
             less reliable estimates or extra assumptions, it is given as "(der)".  Such values are also marked with "(" in
             binast_d.txt.  Where there is no entry for a given parameter (either because it has not been estimated, or because it
             is just an auxiliary parameter needed for derivation of other parameters in cases where they have not been estimated more directly
             and its value was not needed to be specified for the particular object), it is marked with "-tbd-".

The columns have fixed widths, and each is preceded by a space character to assure the columns are separated.
The first line (same in each file) gives the widths of the columns with individual parameters; that facilitates
reading of the files with a program.
