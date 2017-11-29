Paper Notes
===========

Add notes, references, ideas etc for the paper here.

[Asteroid Lightcurve Photometry Database](http://alcdef.org/)
[Minor Planet lightcurve database](http://www.minorplanet.info/lightcurvedatabase.html) (is this just a repackaging of the ALCD?)



Review papers
--------------
* [Asteroids: Recent Advances and New Perspectives](http://adsabs.harvard.edu/abs/2015aste.book....3M) <br>
* [Asteroid Rotations](http://adsabs.harvard.edu/abs/2002aste.book..113P) <br>
* [The Yarkovsky and YORP Effects](http://adsabs.harvard.edu/abs/2015aste.book..509V) <br>
* [Asteroid Systems: Binaries, Triples, and Pairs](http://adsabs.harvard.edu/abs/2015aste.book..355M) <br>
* [Asteroid Models from Multiple Data Sources](http://adsabs.harvard.edu/abs/2015aste.book..183D) <br>
* [Formation and Evolution of Binary Asteroids](http://adsabs.harvard.edu/abs/2015aste.book..375W) <br>

Methods papers
--------------

Add papers describing currently used methods here

* [Waszczak iPTF summer school presentation](http://phares.caltech.edu/iptf/iptf_SummerSchool_2014/slides/waszczak_asteroid_lightcurves.pdf) <br>
This is worth a look - I only skimmed it so far, but I bet it has more references that will be useful too.

* [Warner, Harris, Pravec 2009](http://www.sciencedirect.com/science/article/pii/S0019103509000566?via%3Dihub) The asteroid lightcurve database <br>
The asteroid lightcurve database paper. This describes the contents of the ALCD.

* [On the maximum amplitude of harmonics of an asteroid lightcurve](http://adsabs.harvard.edu/abs/2014Icar..235...55H)  <br>
* [Photoelectric Observations of Asteroids 3, 24, 60, 261, and 863](http://adsabs.harvard.edu/abs/1989Icar...77..171H)  <br>
* [Photometric survey, modelling, and scaling of long-period and low-amplitude asteroids](https://arxiv.org/abs/1711.01893)  <br>
* [Fast and Slow Rotation of Asteroids](http://adsabs.harvard.edu/abs/2000Icar..148...12P)  <br>
* [Modeling of lightcurves of binary asteroids](http://adsabs.harvard.edu/abs/2009Icar..200..531S)  <br>
* [Rotational properties of asteroids, comets andTNOs](http://adsabs.harvard.edu/abs/2006IAUS..229..439H)  <br>


### Lightcurve ###

* [Erasmus et al 2017](http://iopscience.iop.org/article/10.3847/1538-3881/aa88be/meta) Characterization of Near-Earth Asteroids Using KMTNET-SAAO <br>
Short time coverage of NEOs, used L-S to determine periods (and machine-learning for taxonomy). 

* [Vaduvescu et al 2017](https://link.springer.com/article/10.1007%2Fs11038-017-9506-9) The EURONEAR Lightcurve Survey of Near Earth Asteroids <br>
 ... and references therein! (the introduction is useful for more pointers). Uses Canopus (software package by Warner) which uses FALC (Fourier analysis of Lightcurves) algorithm developed by Harris et al. (1989)

* [Harris et al 1989](http://www.sciencedirect.com/science/article/pii/0019103589900158?via%3Dihub) Photoelectric observations of asteroids 3, 24, 60, 261, and 863 <br>
The FALC reference

* [Distribution of shape elongations of main belt asteroids derived from Pan-STARRS1 photometry](http://adsabs.harvard.edu/abs/2017arXiv170905640C) <br>
* [Distribution of spin-axes longitudes and shape elongations of main-belt asteroids](http://adsabs.harvard.edu/abs/2016A%26A...596A..57C) <br>
* [Shape and spin distributions of asteroid populations from brightness variation estimates and large databases](https://arxiv.org/abs/1703.07178)  <br>
* [Spin rate distribution of small asteroids](http://adsabs.harvard.edu/abs/2008Icar..197..497P) <br>
* [Looking a gift horse in the mouth: Evaluation of wide-field asteroid photometric surveys](http://adsabs.harvard.edu/abs/2012Icar..221..226H) <br>
* [ASTEROID ROTATION PERIODS FROM THE PALOMAR TRANSIENT FACTORY SURVEY](https://arxiv.org/abs/1201.1930)  <br>
* [Photometry and Spin Rate Distribution of Small-Sized Main Belt Asteroids](https://arxiv.org/abs/0811.1223)  <br>
* [The asteroid lightcurve database](http://adsabs.harvard.edu/abs/2009Icar..202..134W)  <br>
* [Using sparse photometric data sets for asteroid lightcurve studies](http://adsabs.harvard.edu/abs/2011Icar..216..610W)  <br>
* [ASTEROID LIGHTCURVES FROM THE PALOMAR TRANSIENT FACTORY SURVEY:
ROTATION PERIODS AND PHASE FUNCTIONS FROM SPARSE PHOTOMETRY](http://iopscience.iop.org/article/10.1088/0004-6256/150/3/75/meta)  <br>

### More Detailed Lightcurve / Shape ###

* [Cibulková 2016](https://www.aanda.org/articles/aa/abs/2016/12/aa29192-16/aa29192-16.html)Distribution of spin-axes longitudes and shape elongations of main-belt asteroids <br>
All-sky survey data turned into shapes/spin axes from sparse photometry. Probably not so relevant right now, actually - but maybe for later.

* [Hari Nortunen, Mikko Kaasalainen 2017](https://arxiv.org/abs/1710.06397) LEADER: fast estimates of asteroid shape elongation and spin latitude distributions from scarce photometry <br>
We should definitely take a look at this paper. Finds spin and shape distributions from large databases of observations -- not clear to me from abstract if it's doing this for entire population in statistical fashion or actually finding some values for individual objects (?). Maybe involves sparse lightcurve photometry, rather than dense sampling. 

* [Cibulková et al 2017](https://arxiv.org/abs/1709.05640) Distribution of shape elongations of main belt asteroids derived from Pan-STARRS1 photometry <br>
Looks at a/b and beta distributions of objects from PS1 + others. Uses "LEADER" to determine the distributions of these. 

* [Marciniak et al 2017](https://arxiv.org/abs/1711.01893) Photometric survey, modelling, and scaling of long-period and low-amplitude asteroids <br>
Definitely more in depth modeling than we need for our work - this is going all the way to shapes and spin poles (light curve inversion using SAGE).  Something to think about for the future though, and could be a good reference for gathering light curve / photometry information. 


Download papers
--------------
To download:

https://www.dropbox.com/s/dlyxjlm3u8dnka5/light_curve_references.tar.gz?dl=0

Gaussian Process Papers
-----------------------

Documents describing Gaussian Processes

* A minimal introduction is [here](https://www.robots.ox.ac.uk/~mebden/reports/GPtutorial.pdf)
* [This book](http://www.gaussianprocess.org) is free, but very computer science-y



