FingerprintCalculator
=====================

Collection of python classes for the generation of circular Morgan fingerprints.

The classes can be used independently in python, or
used with the script in which they are integrated.
For that, type: ./FingerprintCalculator.py -help

Fingerprints are calculated in both binary and count format.
Additionally, the fingerprints can be:

1. Hashed: fixed length defided by the user. Collisions might appear.
2. Unhashed (keyed). Each position in the fingerprints corresponds to only one chemical substructure.
The classes as well as the script permit to know which chemical substructures are mapped 
at each position, thus enabling the explanation of the fingerprints and the deconvolution of the 
chemical space if used in predictive modelling studies (e.g. http://www.jcheminf.com/content/7/1/1/abstract).

Please contact me for further details, feedback or help. Thanks !

Example
=====================

./FingerprintCalculator.py -bits 256 -radii 1 2 3 -mols compounds.sdf -output compound_fps -RDkitPath $RDBASE   
