                       ASTEROID LIGHTCURVE DATA BASE (LCDB)
                            Revised 2018 October 7


*******************************************************************************
                              SPECIAL NOTICES
*******************************************************************************
See Section 5 "SPARSE DATA AND WIDE-FIELD SURVEYS"

Mapping in LC_DAT file changed starting with HBAND column to allow
2-character field, e.g., SG instead of g, since upper/lower case
could be confusing.

The format of the LC_BINARY table changed signifiicantly for the Spring 2016
release. See Section 4.3 for more information.

Please review the column mappings for ALL tables. A number were modified for 
the 2016 December release.

Starting with the 2017 annual release, a "basic" summary table is included 
with the frequency-diameter plots based on that summary table. See Section 
2.0.0 for important information on this file and its preferred use.

Releases on/after 2017 August 26 include a new file: LCDB_DAMIT.TXT. See the 
discussion under 2.1.1.1 below.

Releases on/after 2017 November 10 include a new file: LCDB_DIAMETERS_PUB.TXT.
See the discussion under 2.1.0. Also, the LC_COLORINDEX_PUB.TXT file now 
contains all entries, not just those that have a lightcurve period and/or 
amplitude.

Updated Quality (U) Code definitions. See section 4.1.2.

*******************************************************************************


1.0.0  INTRODUCTION
--------------------------------------------------------------------------------

The Asteroid Lightcurve Data Base (LCDB) is a set of files generated from a
dBase IV data base that includes information directly and indirectly obtained
from observations made to determine the period and/or amplitude of asteroid
lightcurves. The information is taken from numerous journals and other sources.

It's main purpose is to provide a central location for basic information about
asteroid rotation rates and related information that can be used in statistical
studies involving a few or many parameters. Some of the data are obtained
directly from the observations while other data are inferred or calculated based
on orbital characteristics, assumed class, etc. Sections below explain in detail
which data are direct and indirectly obtained and, for the latter, their
derivation.

N.B. Even direct data should be confirmed by reference to the original works
     whenever possible. Indirect data are provided for information purposes
     only. They should NOT be used in critical studies.


1.1.0  AUTHOR INFORMATION
-------------------------

These data files are maintained by Brian D. Warner (Palmer Divide Observatory/
MoreData!, Alan Harris (MoreData!), and Petr Pravec (Astronomical Institute, 
Czech Republic).

For basic information on the database or updated versions of the report files,
contact:

  Brian D. Warner
  Center for Solar System Studies / MoreData!
  446 Sycamore Ave.
  Eaton, CO  80615
  brian@MinorPlanetObserver.com

For more theoretical details and discussions, contact

  Alan Harris
  MoreData!        
  4603 Orange Knoll Ave.         
  La Canada, CA 91011  USA
  harrisaw@att.net


1.2.0  DISCLAIMER
-----------------

We have made every attempt to keep the data up to date and correct. However, we
know that there is the possibility for omissions or errors. Please let us know
of any corrections or additions by sending email to one of the below.

  Brian D. Warner   brian@MinorPlanetObserver.com
  Alan Harris       harrisaw@att.net


2.0.0  DATA FILES
--------------------------------------------------------------------------------

The data files in this release consist of two primary sets. The general files
and the Binary Asteroid files provided by Petr Pravec, Astronomical Institute,
Czech Republic.

The general text files contain the most critical information from all the main 
data tables (lookup tables are not included). 

The main exception is the F-D_<Release Date>_PUB_BASIC.TXT file. This contains 
the data used to generate the frequency-diameter plots included in the release. 

THIS FILE DOES NOT INCLUDE U <= 1+, WIDE-FIELD, OR SPARSE-DATA OBJECTS. 

We chose to limit this file to preclude data that we consider less than 
fully-reliable for rotation or other studies.

We STRONGLY urge that the F-D*_PUB_BASIC.TXT file be used as the base for 
initial studies and/or used as a background upon which to plot new results 
to show that they either fit the background, or reveal something profoundly 
new.

Should there be a need for the excluded data, they can be extracted from the 
LC*.TXT files using the column maps and descriptions for each table.


2.0.1  AVAILABILITY OF RAW DATA FILES
-------------------------------------

The raw databases files are not available since they are in a proprietary 
format. We can, on a limited basis, provide CSV files generated from the 
database files. Contact the authors detailing the reaons for needing the 
files and if there are any special requirements to be considered.


2.1.0  GENERAL FILES
--------------------

README.TXT             This file of introductory information.

LC_DAT_PUB.TXT         The lightcurve data file, with data from individual
                       references

LC_SUM_PUB.TXT         Summary data, one line per asteroid, no references.

LC_REF.TXT             Reference list to all published lightcurve data.

LC_AMBIGUOUS_PUB.TXT   Separate summary/detail file for asteroids with ambiguous
                       periods

LC_BINARY_PUB.TXT      Separate summary/detail file for binary asteroids

LC_NPA_PUB.TXT         Separate summary/detail file for "tumbling" asteroids

LC_SPINAXIS_PUB.TXT    Separate summary/detail file for asteroids with spin axis
                       results

LC_NOTESEX_PUB.TXT     Additional notes for Summary/Details entries

LC_COLORINDEX_PUB.TXT  Separate summary/detail file for color indices of
                       asteroids

LC_DIAMETERS_PUB.TXT   Separete summary/detail file for diameters found by 
                       thermal, occultation, radar, or other observations. 


2.1.1  SPIN AXIS CATALOGS
-------------------------

As noted above, the LC_SPINAXIS_PUB.TXT file stores information about the spin
axis properties (ecliptic coordinates and sidereal periods). A more complete and
thorough catalog is maintained by Kryszczynska et al. at the Poznan Observatory 
in Poland.

That catalog can be accessed via

  http://vesta.astro.amu.edu.pl/Science/Asteroids/


Josef Durech also maintains a list of spin axis solutions, with shape models and
data files.  His site is at

   http://astro.troja.mff.cuni.cz/projects/asteroids3D/

It should be noted that the favored DAMIT pole may differ from the one in the 
original reference. This is likely because the original data were used by Durech 
and associates in a new analysis. Their revised result replaced the original 
instead of indicating a new result under a different reference.

Because of the complexities of cross-checking the LCDB vs. DAMIT vs. the original 
result, the LCDB does not directly include any DAMIT results, i.e., there are 
no entries with DAMIT being the author reference. When and if revised results are 
published in one of the journals, those results will included in the LCDB.

2.1.1.1  DAMIT CATALOG (LCDB_DAMIT.TXT)
---------------------------------------
The public release of the LCDB text files now includes a file that summarizes 
the DAMIT catalog by merging the data and reference files from the site. The 
LCDB_DAMIT file is updated with the latest entries from the DAMIT site for 
each public release. In addition, the summarized data are included in MySQL 
tables used on the ALCDEF (alcdef.org) and MinorPlanet.info web sites and 
included in some of the reports generated on those sites.

Field      Header     Format        Pos   Notes
--------------------------------------------------------------------------------
Number     NUMBER     I7            1-7   Blank if no MPC assigned number
Name       NAME       A20          9-28   Name or designation
Citation   CITATION   A20         30-64   Short reference, e.g., Torppa et al. (2003)
L1         L1         I3          66-68   Ecliptice longitude (J2000) of pole
B1         B1         I3          71-73   Ecliptice latitude (J2000) of pole
Period     PERIOD     F12.8       77-88   Sidereal period, hours
Quality    QF         A5          91-95   Quality rating of pole (see DAMIT docs) 
Reference  REFERENCE  A43        97-139   Detailed reference (see NOTES)
Page       PAGE       A10       141-150   Page number(s) in reference (see NOTES)

Notes:
--------------
The DAMIT catalog has two reference ID columns. The primary is RefID; the other is 
RefIDObsolete, which is mainly for an entry that has since be supplanted by a newer 
model. If a Citation entry is followed by '*', this means that there were multiple 
references in the RefID column, with the first entry used for the line in the text.
If the Citation entry is followed by '#', this means that the line is based on the 
first of two or more refernces in the RefIDObsolete column. There is no marker to 
indicate if a single reference for the line used the RefID or RefIDObsolute column.

There is some inconsistency in the DAMIT entries regarding the full reference. In 
particular, whether or not the page number was part of the value in the Publication 
column. For these reason, the values in the Publication and Page columns were taken 
at face value, with the exception that journal names were made consistent. 

For example, there were several variations on "Astronomy & Astrophysics." These 
were all changed to "Astron. Astrophys."


2.2.0  BINARY ASTEROID FILES
----------------------------

The following files are maintained by Petr Pravec. These files are considerably
more detailed than what's in LC_BINARY_PUB and so are better suited for critical
studies.

BINARY_README.TXT  Separate README pertaining the binary asteroid
BINASTD_PUB.TXT    The best estimates of compiled parameters
BINASTE_PUB.TXT    Uncertainties of the estimates in BINASTD_PUB.TXT
BINASTM_PUB.TXT    References and notes for the compiled estimates
BINASTR_PUB.TXT    Information on each of the estimates, e.g., their derivation

These files are available at

  http://www.asu.cas.cz/~asteroid/binastdata.htm



2.2.1  ABOUT BINARY DATA
------------------------

The data in the LC_BINARY_PUB.TXT file are by no means exhaustive. They are
meant to provide a quick overview of the parameters for a binary (or multiple)
system.

The reader is urged to consult the original journal articles for more complete
details.

A good summary page with links to those journal references is the web site run
by by Wm. Robert Johnston

  http://www.johnstonsarchive.net/astro/asteroidmoons.html

Another summary database, the Virtual Observatory Binary Asteroids Database (VOBAD),
is maintained at UC Berkeley/SETI Institute by Franck Marchis, Brent Macomber, et al.

  http://cilaos.berkeley.edu/PHP_scripts/VOBAD/VOBAD_portal.html



3.0.0  LCDB DATA
--------------------------------------------------------------------------------

The original lightcurve database was a simple text file with a structure that
tried to minimize disk space requirements. That served well for many years.
However, the rapidly growing number of lightcurves being reported as well as the
file's format not being able to accommodate some data prompted a change starting
in mid-2006.

The foremost change was converting to a relational database that includes
numerous tables and has SQL search capabilities. This allows for not only easier
maintenance of the database but for generating reports in a new way that are
more informative, complete, and constant in formatting.

The format of the report files has changed from previous years. Some significant
changes were made for the Spring 2016 release to include additional output fields.
Third-party parsing programs will need to be modified to handle the changes. See
the revised column mappings in Section 4.

The sections below provide the formatting for each field and the reports. This
should make rewriting old parsing programs, or new ones, fairly straightforward.


3.1.0  DIRECT DATA
------------------

Data that is obtained directly from photometric observations includes

  1. Rotation period (usually synodic).

  2. Amplitude.

  3. Absolute magnitude, H, and phase slope parameter, G, when determined by
     using reduced magnitude versus phase angle data.

  4. Binarity due to mutual events, i.e., occultations and eclipses. In such
     cases, the rotation period of the primary and orbital period of the
     satellite and the amplitude of the primary lightcurve are the usual direct
     results. The size ratio can be computed from the depth of the events.
     For more details on binary lightcurve analysis, see Pravec et al., 2006,
     Icarus 181, 63-93.

  5. Color indices.

  6. Diameter if based on stellar occultation or adaptive optics/radar. Radar
     diameters can also be considered indirect, depending on how the diameter
     was determined.

  7. Taxonomic class.



3.1.1  SYNODIC VERSUS SIDEREAL PERIOD
-------------------------------------

The synodic period depends on viewing aspect and the rate of motion of the
asteroid across the sky.  An expression for the magnitude of the expected
difference between the sidereal period and synodic period based on the Phase
Angle Bisector (PAB) is

   DeltaP = [d(PAB)/dt] * P^2

Where  DeltaP       difference between synodic and sidereal periods, in units
                    of the rotation period (usually hours).

       [d(PAB)/dt]  angular rate of change in the Phase Angle Bisector in
                    inverse units of the rotation period, e.g., in units
                    of cycles/hour

       P            the rotation period of the asteroid in the same units of 
                    time, e.g., hours/cycle

For example, assume an asteroid with a rotation period of 8 hours observed when 
the PAB is changing by 0.05 deg/day (typical for an main-belt asteroid at 
opposition), then the sidereal-synodic difference is

  DeltaP = [0.05 / 360.0 / 24.0]  *     (8 ^ 2)       =   0.00037 hr.
             cycles / hour          hours^2/cycles^2     hours / cycle

The synodic-sidereal difference can be either positive or negative, and can
exceed the value given by this expression for near-polar aspects, but the
expression gives a reasonable estimate of the magnitude of the expected
difference.

In most cases, the period given in the summary and details lines are synodic and
not sidereal. An 'S' flag (see notes below) indicates that the period is
sidereal. There are many entries in both tables that do not carry the 'S' flag
when they should. This is part of the legacy nature of the data after converting
the files to the new data base, i.e., the old format did not allow for
indicating the period was one type or another. We hope to update these and other
legacy values that now have qualifying flags in future versions.

For most studies, the difference between sidereal and synodic period is not
significant.


3.1.2  INDIRECT DATA
--------------------

Indirect data is that obtained by calculation and/or assumption.

N.B. THESE DATA ARE PROVIDED FOR CONVENIENCE ONLY, AND SHOULD NOT BE USED FOR
     STATISTICAL STUDIES, NOR ASSUMED TO BE ACCURATE AND CURRENT.


  1. Diameter, H, and albedo (Pv, Pr, etc.).

     The relationship between H, diameter, and albedo is:

         D = (1329km) * 10.0^(-0.2*H) / sqrt(albedo)
     or

         logD (km) = 3.1235 - 0.2H - (0.5 * log(albedo)).

     (See Pravec and Harris, 2007. Icarus 190, 250-259).

     The value of H is usually known, though not always accurately, based on
     photometric observations. If direct data are available for D and/or Pv,
     then the above relationships can be used to derive a missing quantity.

     Data from the SIMPS study (Tedesco, E.F., Noah, P.V., Price, S.D.: 2004,
     IRAS Minor Planet Survey. IRAS-A-FPA-3-RDR-IMPS-V6.0. NASA Planetary Data
     System) is used when available and no overriding data are available. If a
     newer value of H than that used by SIMPS is available, the diameter and
     albedo are re-computed based on Harris and Harris (1997, Icarus 126,
     450-454).

     If the diameter was based on an assumed albedo and H is revised, the albedo
     is held constant and the diameter re-computed using the above formulae.

     If the diameter was determined by radar, resolved imagery, etc. and a new
     H is available, the diameter is held constant and the albedo is
     re-computed.


  2. Taxonomic Class, orbital class, and albedo

     These three values can have a complex relationship when the class and
     albedo are not directly obtained.

     When spectroscopic or other data are available that can determine the
     taxonomic class exactly but no direct albedo data are available, the albedo
     can be assumed. This albedo can, in turn, be used to find the estimated
     diameter and/or H if those quantities are not directly known.

     Flags in the Summary table indicate the source of these fields, including
     if they are assumed based on a combination of available data. The table
     below shows the assumed values used barring any direct data.

     The family/group is based on orbital parameters. These are mostly for
     informational purposes only since the definition for some groups or
     families is "fuzzy", at best. Note the distinction between a family and
     group. A family is a set of asteroids with a common parent body. A group
     is a set of asteroids with common orbital characteristics. Members of a
     family will almost always be in the same group, but members of a group may
     not necessarily be of the same family.

     For a more complete discussion of this topic, see Warner et al. (2009)
     Icarus 202, 134-146. This paper includes new default albedo and G values 
     for given taxonomic classes based on the data in the LCDB.


     Group/Family                Source/orbital parameters                      Class  Pv   Note
     -------------------------------------------------------------------------------------------
     Baptistina                  Bottke el al                                     C   0.057
     Centaur                     5.5 < a <= 30                                    C   0.057
     Centaur Comet               Centaur w/comet behvaior                         C   0.057   9
     Comet-like orbit            Q > 5.0                                          C   0.057   4
     Comet                       exhibits coma and/or tail                        C   0.057
     Eos                         2.99 < a < 3.03, 0.01 < e < 0.13,  8 < i  < 12   S   0.14
     Erigone                     2.32 < a < 2.40, 0.15 < e < 0.22,  4 < i  <  6   C   0.057
     Eunomia                     2.53 < a < 2.72, 0.08 < e < 0.22, 11 < i  < 16   S   0.21
     Flora                       2.15 < a < 2.35, 0.03 < e < 0.23, 1.5 < i < 8    S   0.24
     Hilda                       3.7  < a < 4.2,  0.07 < e,              i < 20   C   0.057             
     Hungaria                    1.78 < a < 2.0, e < 0.18, 16 < i < 34            ES  0.3     1
     Karin                       Nesvorny, private communications                 S   0.26
     Koronis                     2.83 < 2.91, e < 0.11, i <= 3.5                  S   0.24
     Main belt - inner           a < 2.6                                          S   0.20    5
     Main belt - inner comet     a < 2.6                                          C   0.057   9
     Main belt - middle          2.6 < a < 2.7                                    SC  0.10
     Main belt - middle comet    2.6 < a < 2.7                                    SC  0.10    9
     Main belt - outer           2.7 < a < 5.0                                    C   0.057
     Main belt - outer  comet    2.7 < a < 5.0                                    C   0.057   9
     Mars crosser                1.3 < q < 1.668, Q < 5.0                         S   0.20    5
     NEA (Apollo/Aten/Amor)      q < 1.3                                          S   0.20    5
     NEA (comet)                 q < 1.3                                          C   0.057   9
     Nysa                        2.4 < a < 2.5, 0.12 < e < 0.21, 1.3 < i < 4.3    S   0.20    6
     Phocaea                     2.25 < a < 2.5, e >= 0.1, 18 < i < 32            S   0.23
     Planet Satellite                                                                         7
     Themis                      3.08 < a < 3.24, 0.09 < e < 0.25, i <= 3         C/B 0.08
     TNO/KBO                     a > 30                                           C   0.1
     Trojan - Jupiter            5.05 < a < 5.4                                   C   0.057
     Trojan - Mars               similar to Mars                                  S   0.20    5
     Trojan - Neptune            similar to Neptune                               C   0.057
     Trojan - Saturn             similar to Saturn                                C   0.057   2
     Trojan - Uranus             similar to Uranus                                C   0.057   2
     Vestoid                     2.26 < a < 2.48, 0.03 < e < 0.16,  5 < i < 8.3   S   0.20    3

     NOTES
     -----
     1  Pv = 0.3 is a compromise value when no taxonomic information is
        available, since the Hungarias are both a family (common parent, E/X
        class, pv = 0.4) and group (similar orbits, S class, Pv = 0.20).

     2  None known and not likely due to perturbations by giant planets, both
        interior and exterior.

     3  Higher albedo (0.4) assigned only if in V class determined (SMASS, etc).
        Otherwise, class is 'S' and Pv = 0.20 on presumption that the object
        shares orbital characteristics but not parent body.

     4  Barring any other classification that meets Q > 5, the orbit is
        classified as "comet-like."

     5  The default Pv = 0.20 ± 0.07 for S-type objects was derived from the
        geometric mean of all S-type objects in the LCDB with known albedos
        (usually SIMPS).

     6	The Nysa orbital space is polluted by a large portion of objects that
        are not true members of the Nsya-Hertha family and even then, the true
        members are heterogeneous in nature. For this reason, we elected to
        treat the Nysa space the same as other unassociated inner main belt
        objects by using a default of S class and Pv = 0.20.

     7  The Planet satellite group includes small and distant natural satellites of the major
        planets, e.g., Himalia Jupiter VI. They are included since a number of
        these smaller bodies resemble asteroids in size and rotational properties.

     9  Several orbital groups have been given a "comet" subclass. This is for
        objects within the given orbital class that have shown cometary activity,
        e.g., so-called "main belt comets." The Comet class is reserved for "true"
        comets, e.g., P/Encke and the Comet-like orbit class is still reserved for
        objects that have very elongated orbits but have shown no signs of cometary
        activity. These "comet" subclasses are given default taxonomic class of 'C' 
        and albedo of 0.057.


     Bottke, W.F., Vokrouhlicky, D., Nesvorny, D. (2007). Nature 449, 48-53, and
     private communications.


  3. Color Index applied to H

     Color index is not generally assumed or entered into the LCDB. However,
     some times the value of H was found in a photometric band other than V,
     e.g., Cousins R. In that case, and if the value is used to override the H 
     given by SIMPS or the MPCORB file (Minor Planet Center) in the Summary
     line, H is transformed to the V band. When the color index is not directly
     available, these values are used to transform the measured H value

       V-R  0.45
       B-V  0.80
       V-r' 0.23

     When such a transform is used in the Summary line, whether or not based on
     an assumed value, the H value has the 'T' (transformed) flag.



4.0.0  LCDB FILE DESCRIPTIONS
--------------------------------------------------------------------------------

The following sections describe the specific files that are part of the LCDB
release, not including the Binary files from Petr Pravec. A column map is
provided for each report that is part of the release. In addition, sub-sections
describe the meaning of the flags that qualify various fields in the report.


4.0.1.  FILE STANDARDS
----------------------

All files are simple ASCII text with Windows (PC) line terminators (CRLF).

Unless stated otherwise, all files are single-space delimited, i.e., there is a
one-character blank column between the end of one field and the start of the
next.



4.1.0  LC_SUM_PUB AND LC_DAT_PUB
--------------------------------

These are the primary files in the LCDB release. Both have the same column
mapping. They show both direct and indirect data, the most important being the
direct data of lightcurve period and amplitude along with our assessment of the
quality of the period solution. The latter is expressed by the U code, which is
described in detail below.

In the LC_DAT_PUB file, a Summary line is followed by one or more Details lines.
The Summary line is our best determination of the primary information for a
given object based on the data in the Details lines. For example, where several
periods are available, the Summary line uses the one that we consider the most
likely to be correct. Sometimes that value may be an average of the available
values.

The Details lines differ in formatting only in that the number and name of the
asteroid are replaced by a "short reference" name that will be found in the
LC_REF_PUB references file.


4.1.1  MULTIPLE DETAILS TABLE ENTRIES
-------------------------------------

In some cases, there is more than one detail line under a given asteroid with
the same publication reference. This is deliberate in order to allow statistical
studies of lightcurve amplitude versus phase versus class (albedo).

For example, if a single publication reports the lightcurve behavior for an
asteroid where the synodic period and/or amplitude of the curve changed
significantly during the course of the observations, the Details table includes
an appropriate number of entries. Those entries will "split out" the results so
that the period and/or amplitude can be tied to a specific (though maybe only
approximate) set of PAB or Phase values. A good example would be a paper
reporting observations of an NEA asteroid over several weeks where the amplitude
of the curve went from 1.1 to 0.3 magnitudes over the range of observations.

In most cases, splitting the results into distinct sets was not difficult, e.g.,
the asteroid was observed on one night at one-week intervals. In some cases, the
split was not so distinct. In this case, compromises were made in order to avoid
having an excess of multiple entries while still retaining sufficient resolution
of the variations versus time.

A variation on the above is if the author(s) forced the data from several blocks
of dates to fit a fixed period solution. Here, the period will be the same for
all entries, though the amplitude may change. In this case, the period is left
blank. The U code is assigned for each lightcurve based on the presumption that 
the fixed period is correct, i.e., it is based on the quality of the fit of the 
data to the presumed period. The main point of interest is the amplitude for the 
reasons given above. Other information that was derived based on the given block 
of data, e.g., a value for H, G, or a color index, will be included within that
Details record as well so that it's clear which block of data was used to derive
the given values.



4.1.2  U (QUALITY) CODE
-----------------------
The U code provides our assessment of the quality of the period solution, not
necessarily of the data per se. The uniqueness of the solution, while an 
important factor, is not the sole consideration in making an assessment. 

Depending on the specifics for a given asteroid, a good period solution can be 
obtained by using a large amount of lesser quality data about as well as using 
less data that is of higher quality. Many factors come into play making the 
assessment. The table below gives the general outline of the criteria used,
going from highest to lowest rating.


  3   The lightcurve is completely unambiguous in terms of period, i.e., 
      there are no cycle ambiguities or possible solutions with single, 
      triple, or other number of extrema. The coverage of the entire rotation 
      phase is to the degree than any remaining small gaps can be confidently 
      interpolated.

  3-  A unique period determination, but possibly some moderate gaps in 
      coverage, enough so that interpolation of the entire curve is not 
      certain, but not enough to allow any other solution. Also, while 
      the period is "unique", it may still be wrong by a small amount due 
      to outliers, large scatter in relation to the lightcurve amplitude, 
      lightcurve evolution over the range of observations or tumbling that 
      leads to imperfect overlapping at some parts of the lightcurve, 
      and/or not finding the period to a high enough precision.

  2+  It is unlikely but not impossible that the period is in error due to 
      cycle counts or alternate numbers of extrema per cycle, and no more 
      than moderate gaps in coverage (as in U = 3-).

  2   Result based on less than full coverage, so that the period may be wrong 
      by 30 percent or so. Also, a quality of 2 is used to note results where 
      an ambiguity exists as to the number of extrema per cycle or the number 
      of elapsed cycles between lightcurves. Hence the result may be wrong by 
      an integer ratio.

  2-  Period and total amplitude not firmly established. For example, a single 
      night coverage of about half a cycle including a maximum and minimum, 
      but not enough to actually derive an accurate period.  This is the 
      minimum reliability code that we accept for statistical analysis. 

  1+  Similar to U = 2-, but with less amplitude so that it is not absolutely 
      certain that the variations are true rotational variation and not due to 
      noise, etc.

  1   May be completely wrong. What is interpreted as rotational variation 
      may be just noise, calibration error, etc.

  1-  Probably wrong. A lightcurve that may be completely wrong (as in U = 1) 
      but, in addition, the claimed period is very unlikely, e.g., a large 
      object with a claim of P < 2h. 

  0   Result later proven incorrect. This appears only on records of 
      individual observations.

      It is important to keep in mind that U = 0 does not necessarily mean
      that the data for a given lightcurve are of low quality. The only 
      interpretation that should be inferred is that the -reported solution- 
      has been determined, perhaps from subsequent data, to be incorrect so 
      that not even the loose constraints of U = 1 or U = 2 can be used. For 
      the most part, U = 0 will be used very sparingly and the previous U 
      rating (unless 3) will be retained to avoid the false impression that 
      the data are of limited or no use. 


  N.B. Until the intermediate release in 2008 November, the LCDB also used a
       value of '4' for the U code, which indicated that a pole solution had
       been reported. This is no longer the case since, in the past, there have
       been cases where a 4 was assigned because there was a pole solution given
       but the best available period solution was no better than 2.

       The period solution quality is now indepedent of any pole solution. A
       separate "Pole" flag in the Summary and Details reports is used to
       indicate that a pole solution has been reported. The LC_SPINAXIS_PUB file
       includes more details and its own quality code assignment.


Assignment of the refined ratings using a + or - is a work in progress as we 
catch up with almost 20 years of data entry. Therefore, not all U code ratings 
will match what we would give under current rules and are subject to change.

If a summary line contains U = 3 it should be understood that, by implication, 
any detail lines with U < 3 are superceded in terms of period determination.  
We will generally retain the U rating assigned to a prior U < 3 result if the 
period derived lies close within its undertainty range to the later adopted 
U = 3 period, as an indicator of the quality of the observations reported.  

In the case where U < 3 in the summary line, then one or more detail lines may
be assigned U = 0, but only if it is firmly established that, while we
cannot say for certain what is the correct period, we can say with certainty
that the period claimed in the reference is wrong (i.e., discordant with 
subsequent data). 

Some Details lines, and even some Summary, may not contain a U code rating. This
is deliberate and can be for several reasons.

  1. The available data do not include a lightcurve, therefore, it is not
     possible to give a rating to the curve. In some cases, where the results
     are reported by observers whose standard of work is known to be of
     sufficient quality, we may assign an interim U code, usually 2, until a
     lightcurve or the data are available.

  2. In the case where several results are published for a given object in the
     same reference, we will assign a U code rating for the "best" available
     data and include only new information for that given Details record, e.g.,
     see section 4.1.1, "Multiple Details Table Entries."

There are instances when data are available they do not reasonably define a 
period or even constrain a range in which the period lies. In addition, the data 
may not be able to provide any reasonable indication of the amplitude. In these 
cases, the entry in the Details record will have the reference to the work and 
blank U code rating and no period or amplitude. 

The Summary line may also have no period and/or amplitude as well as U = <blank>.
This occurs when none of the Detail lines, even if they have some or all of the
information, is deemed insufficiently reliable so as to put the information in the
Summary line. This is done to show that there are data available for the object 
but that they may be of very limited use.




4.1.3  LC_SUM_PUB AND LC_DAT_PUB COLUMN MAPS
--------------------------------------------

In the LC_DAT_PUB file, there is a blank line between the last detail line for a
given object and the Summary line for the following object.


Field      Header   Format        Pos   Notes
--------------------------------------------------------------------------------
Number     NUMBER   I7            1-7   Blank if no MPC assigned number
EntryFlag           A1              9   * new record since last pubic release (see notes)
Name       NAME     A30         11-40   Summary: Name or designation
                                        Details: Publication reference
Desig      DESIG    A20         42-61   Summary: MPC Designation (see notes)
                                        Details: Appxoimate mid-date of observations
Family     FAM      A8          63-70   Family/Group association
CSource             A1             72   Source of taxonomic class
Class      CLASS    A5          74-78   Taxonomic class
DSource             A1             80   Source of diameter
DFlag               A1             82   Diameter flag
Diameter   DIA.     F9.4        84-92   in km
HSource             A1             94   Source of H value
H          H        F5.2       96-100   H (absolute magnitude)
HBand               A2            102   BVRI,S(UGRIZ), V if Blank
ASource             A1            105   Source of albedo value
AlbedoFlag          A1            107   Albedo qualifier (> or <)
Albedo     A        F6.4      109-114   Albedo
PFlag               A1            116   Period qualifier
Period     PERIOD   F13.8     118-130   Period, in hours
PDescrip   P DESC   A15       132-146   Period description if not numeric value,
                                        e.g., long
AmpFlag             A1            148   Amplitude flag
AmpMin     AMIN     F4.2      150-153   Minimum amplitude of a range.   (See Notes)
AmpMax     AMAX     F4.2      155-158   Maximum amplitude of a range OR (See Notes)
                                        amplitude if no range.
U          U        A2        160-161   Lightcurve Quality
Notes      NOTES    A5        163-167   Qualifying flags for lightcurve record
Binary     BIN      A3        169-171   Binary Flag
Private    PRI      A3        173-175   Y = Unpublished, contact named observer
                                        to request details
Pole       SAM      A3        177-179   Spin Axis and/or Shape Model available

SparseData SD	    A2        181-182   Y = Result based on sparse data (see Sec. 5)

WideField  WF       A2        184-185   Y = Result based on wide field data (see Sec. 5)

Notes:
--------------
Entry flags (column 9)
The once-a-year "official" release of the LCDB will have only of two values in
column 9

  <blank>  The entry was last entered/updated more than one annual release
           prior to the current version.

     *     The entry is new since the last annual release.

     O     (Capital O) The summary record is an "orphan", i.e., no lightcurve data 
           have been reported. However, other data such as color indices, H-G values, 
           etc. may have been reported and so a record is created in the summary 
           table with the hope that lightcurve data will be reported. 

           The only time orphaned summary records are included is in the color index 
           and diameters tables. Otherwise, fewer than 30% of the entries in the 
           color index table and an even smaller percentage of the diameter table 
           entries would be available in the public release.

Intermediate releases, those not fully vetted and so offered "as-is", might
contain other values in column 9, e.g., '@' and '!'. These are "bookkeeping"
marks used during the vetting process. The only thing to be inferred from them
is that the entry has yet to be validated, is subject to change, and may
eventually appear with the '*' flag.

Summary records marked private include only the U code and private flag. The
accompanying Details line(s) include only the author, and the U code. This
allows the LCDB user to know that there is a lightcurve available and its
source. He should contact the author to inquire about the results and
availability of data.

The DESIG field in the Summary line is reserved for the official MPC designation
only. It is not used to record survey/catalog specific IDs. It can be the same
as in the NAME field until the object is given an official name by the IAU.

Asteroids 1 - 332 do not have MPC designations.

The DESIG field in the Detail line gives the mid-date of observations, if available,
in yyyy-mm-dd format, e.g., February 20, 2001 = 2001-02-20. It is this date for
which the phase angle and phase angle bisector values are computed.

For example,

  70030 * Margaretmiller                 1999 CZ1
          Warner 2005f                   2003-10-06
          Warner 2010k                   2003-10-06
          Warner 2010k                   2010-02-15

As opportunities allow, the Details DATEOBS field will be updated to give the
mid-date of observations for a specific data set.

The Min/Max Amplitude values in the Summary table are based only on detail 
lines that have a U >= 2- ratings. If the U code is empty or U <= 1+, the 
detail line min/max amplitudes are not considered.


4.1.4  FIELD (FLAG) CODES USED IN SUMMARY AND DETAIL LINES
----------------------------------------------------------

The flags appear in the data field immediately before the value they qualify. In
most cases, they are a single character.

ALBEDOFLAG (Albedo Qualifier)
-----------------------------
Blank   NONE
<       Less than
>       Greater than


AMPFLAG (Amplitude Flag)
------------------------
Blank   NONE
<       Less than
>       Greater than


ASOURCE (Albedo, SUMMARY)
-------------------------
Blank   Default = measured
A       Assumed (based on class - known or assumed)
C	Caculated from H and Diameter
D       Derived (usually from IR measurements were diameter is measured directly)
M       Measured (radiometric or polarimetric observations)


ASOURCE (Albedo, DETAILS)
-------------------------
Blank   Unknown source
A       Assumed (based on class - known or assumed)
D       Derived (based on H and/or D being known; Harris & Harris method)
L       Taken from Details List entry
S       SIMPS (IRAS Minor Planet Survey)


BIN (Binary Flag)
-----------------
Blank   NONE
?       Possible, usually due to lacking mutual eclipse/occultation events
B       Binary system
M	Multiple system, more than two


If the '?' flag appears in a detail line, the author made a claim that the 
object is or may be a binary. If the '?' also appears in the summary line, 
we regard the claim as credible.

If '?' appears only in the detail line and not in the summary line, we regard 
the binary claim as spurious or unfounded. 


CSOURCE (Class Source)
----------------------
Blank   NONE
A       Assumed (see discussion 3.1.2, "Indirect Data");
B       SDSS-Bus (derived by Bus et al. from SDSS colors)
L       Details List
M       SDSS-MFB (derived by Masi et al. from SDSS colors)
S       SMASS I/II
T       Tholen et al, Asteroids II



DFLAG (Diameter Flag)
---------------------
Blank   NONE
<       Less than
>       Greater than


DSOURCE (Diameter Source)
-------------------------
Blank   NONE
C       Calculated A/H   Based on assumed albedo and MPCORB H - see discussion
                         below
D       Derived A/H      Recalculated from determined albedo and/or H - see
                         discussion below
L       Details List     Taken directly from a details list value
S       SIMPS            Supplementatl IRAS Minor Planet Survey


EntryFlag
---------
Blank   Pre-existing record.
*       New or updated record since the last release.



FAMILY (ORBITAL GROUP)
----------------------------------
Blank   No group or family determined
BAP     Baptistina
CEN     Centaur
CEN-C   Centaur - comet behavior (coma and/or tail)
CLO     Comet-like orbit
CMT     Comet
EOS     Eos
ERI     Erigone
EUN     Eunomia
FLOR    Flora
H       Hungaria
HIL     Hilda
KAR     Karin
KOR     Koronis
MB-I    Main belt - inner
MB-IC   Main belt - inner, comet behavior (coma and/or tail)
MB-M    Main belt - middle
MB-MC   Main belt - middle, comet behavior (coma and/or tail)
MB-O    Main belt - outer
MB-OC   Main belt - outer, comet behavior (coma and/or tail)
MC      Mars crosser
NEA     NEA (Apollo/Aten/Amor)
NEA-C   NEA (Apollo/Aten/Amor), comet behavior (coma and/or tail)
NYSA    Nysa
PHO     Phocaea
THM     Themis
TNO     TNO/KBO
TR-J    Trojan - Jupiter
TR-M    Trojan - Mars
TR-N    Trojan - Neptune
TR-S    Trojan - Saturn
TR-U    Trojan - Uranus
V       Vestoid

See the discussion in 3.1.2, "Indirect Data".


GSOURCE (G = Phase Slope Parameter, SUMMARY)
--------------------------------------------
Blank   None  (assumed to be 0.15 in MPCORB)
D       Default (assigned based the assumed value for a given taxonomic class)
L       Taken from Details list entry
G	G12 system
M       MPCORB (usually 0.15)


GSOURCE (G = Phase Slope Parameter, DETAILS)
--------------------------------------------
Blank   None  (default = measured)
A       Assumed (based on class, orbital group, etc.)
C	Calculated
D	Derived
G	G12 system
M	Measured

HBAND (H = absolute magnitude, standard filter)
-----------------------------------------------
Blank   Assumed to be Johnson V
B       Johnson Blue filter
V       Johnson Visual filter.
R       Red filter (presumed to be Cousins R, not Johnson R)
I       IR filter (presumed to be Cousins I, not Johnson I)
SG      SDSS g
SR      SDSS r
SI      SDSS i


HSOURCE (H = absolute magnitude, SUMMARY)
-----------------------------------------
Blank   None or unknown.
A       ASTORB file maintained by E. Bowell, Lowell Obs., Flagstaff
D       Derived from reported Diameter and or albedo
E	Estimated, e.g., derived from a "Bowell orbit" and assumed albedo
L       Taken from Details list entry
M       MPCORB file maintained by the Minor Planet Center
S       SIMPS (IRAS Survey)
T       Transformed


The T (Transformed) flag is used only in those cases where the adopted H value
is based on a details line but the reported value in that line was not in the V
band. If the actual color index is available, it is used to convert the Details
H value to Johnson V. If no color index is availble, default color indices are
used to convert H_R and H_B to H_V:

   V-R  0.45
   B-V  0.80
   V-r' 0.23

These are slightly red of solar color, reflecting the tendency of asteroids in
that direction.


HSOURCE (H = absolute magnitude, DETAILS)
------------------------------------------
Blank   Default = Measured
A       Assumed
C       Calculated from pv and Diameter
D       Derived (usually from IR observations where diameter was measured directly)
M       Measured (direct measurement from reduced mags over range of phase angles)



NOTES (single letter flag(s))
-----------------------------
Blank   NONE
?       Usually tied with 'T' or 'A' flags to indicate uncertainty
-       Tied with T flag. See notes below.
<X>     X is an integer, e.g., 3 for a trimodal curve (three min/max per rotation)
          or 4 quadrimodal curve (four min/max per rotation), etc.
A       Ambiguous period (see LC_AMBIGUOUS_PUB.TXT for details)
D       Period determined by us that differs from that given in the original
          publication
E       Occultation observation (usually when reporting a diameter)
H       Space telescope observations (optical)
I       IR/Thermal observations (e.g., Spitzer)
M       Polarimetric observation
N       No lightcurve published
O       Adaptive optics observation
P       Photographic photometry
R       Radar observation
S	Spectroscopic
T       Tumbling (NPA rotation - see LC_NPA_PUB.TXT for details and notes below)
V       Visual photometry


The 'A' and 'T' notes flags are used to call the reader's attention to the
LC_AMBIGUOUS or LC_NPA reports, repsectively. They should not be taken as stand-
alone information. Instead, consider them footnote numbers in the body of a main
text. The other reports (and original references) are the actual footnotes.

The A flag does not appear in the Summary line unless the Summary line value
itself represents an ambiguous solution, i.e., just because a Details line may
report an ambiguous period does not mean that the Summary period is also
ambiguous.

The T flag currently has four possible qualifiers:

Blank   The asteroid has a PAR < -1, i.e., it is definitely tumbling.
        Example: T

?       Possible tumbler. There is some evidence that the asteroid might
        be a tumbler. It may carry a PAR = 0 to -1. See the discussion 
        for the LC_NPA_PUB table for the meaning of the PAR codes.
        Example: T?

0       The tumbling damping time scale (see Pravec et al. 2014, Icarus) is 
        long enough that tumbling might be expected, but observations are not
        sufficient to substantiate either tumbling or not tumbling, PAR = 0.
        Example: T0

-       The tumbling damping time scale is long enough that tumbling might
        be expected, but observations indicate that the object is NOT tumbling,
        i.e., PAR >= 1.
        Example: T-

+       The tumbling damping time scale is short enough that tumbling would not
        seem likely, however observations indicate that it may be tumbling or
	actually is tumbling.  
        PAR = < 0.
        Example: T+

We include the expanded tumbling notes to call attention to what we consider to
be an important aspect in the study of YORP spin up/down theories. This is done
by making known any asteroids that are or are strongly believed to be tumbling
as well as those that should be and aren't or are and shouldn't be.

The W flag is included so that those doing statistical studies can include or 
exclude the results from these surveys. Such surveys can introduce significant 
biases by "cherry picking" the best results from a large pool and so skew 
overall rotational statistics. See the paper by Warner and Harris (2011, Icarus).


PFLAG (Period Flag)
-------------------
Blank   NONE
<       Less than
>       Greater than
D       No numerical value, see P DESC field description
S       Sideral period, default is no flag and synodic period
U       Uncertain, not the same as ambiguous where one or additional periods are
        reported. For example, the data did not allow finding a definite period
        and so the author(s) reported a "best guess."


4.1.5  DATA SUITABLE FOR ROTATION RATE STUDIES
----------------------------------------------
As noted in Warner et al. (2009, Icarus 202, 134-146), only those objects
with a U code of 2- or greater in the LC_SUMMARY file, i.e., U = 2-, 2, 2+, 3-, 
or 3, should be used for rotational rate studies and, unless there is a 
specific reason otherwise, the summary line period should be used instead of
one of the periods in the details table. 


4.2.0  LC_AMBIGUOUS_PUB.TXT  (AMBIGUOUS PERIODS)
--------------------------------------------------------------------------------

This file includes any record where the Notes flag for a Summary and/or Detail
line indicates an ambiguous period.

There is not always a direct cross-connection between the Summary and Details
entries. For example, it's possible to have a Summary line  without the
ambiguous period flag but one or more of the Details lines to have the flag. In
this case, we judge that the ambiguity has been resolved by subsequent
observations, but retain the ambiguous flag in the detail line for historical
accuracy.

In turn, if the Summary line is flagged as ambiguous, this does not mean that
any of the Details lines are also flagged as such. In that case, it means that
no one solution sufficiently stands out and so the one that is reported on the
Summary line is considered to be only the most probable solution.

The first line for a given object is the Summary line, which contains the number
and name of the object and the adopted period and amplitude. As noted above, the
Details lines(s) may not agree with the Summary line.


4.2.1  LC_AMBIGUOUS_PUB COLUMN MAPPING
--------------------------------------

Field       Header     Format      Pos   Notes
--------------------------------------------------------------------------------
Number      Number     I7          1-7   Summary: MPC assigned number. Blank if
                                                  none.
                                         Details: blank.
EntryFlag              A1            9   * = New entry since last release.
Name        Asteroid   A30       11-40   Summary: Name of object.
ShortRef    Reference  A30       11-40   Details: Publication reference.
Notes       Notes      A5        42-46   Qualifying flags for summary and detail
                                         lines. See the NOTES description in
                                         section 4.1.4.
Period      Per1       F13.8     48-60   Summary: Adopted period.
                                         Details: Most likely period.
AmpMax      Amp1       F4.2      62-65   Summary: Adopted maximum amplitude.
                                         Details: Amplitude for most likely
                                                  period.
Period2     Per2       F13.8     67-79   Details: Second possible period.
Amp2        Amp2       F4.2      81-84   Details: Amplitude for second period.
Period3     Per3       F13.8     86-98   Details: Third possible period.
Amp3        Amp3       F4.2    100-103   Details: Amplitude for third period.
Period4     Per4       F13.8   105-117   Details: Fourth possible period.
Amp4        Amp4       F4.2    119-122   Details: Amplitude for fourth period.
Period5     Per5       F13.8   124-136   Details: Fifth possible period.
Amp5        Amp5       F4.2    138-141   Details: Amplitude for fifth period.



4.3.0  LC_BINARY_PUB.TXT  (BINARY ASTEROIDS)
--------------------------------------------------------------------------------

This file includes those asteroids that are known or suspected binaries. This is
not meant to be a comprehensive compilation of data for binary asteroids. Visit
the URL given in section 2.2.0 for a page that provides more details as well as
links to the original journal articles.

The Summary line gives the number (if any) and name of the object as well as the
adopted period and amplitude.

The Details table lists core data, including the type of binary. There are three
broad categories:

  Fully-asynchronous   (Binary Type Flag: A)  Example: 1509 Esclangona
  The satellite's rotation period is different from its orbital period. In this
  case, the orbital period is given along with the independent rotation period
  and lightcurve amplitude of the satellite, if available.

  Singly-asynchronous  (Binary Type Flag: S)  Example: 5905 Johnson
  The satellite's rotation period and orbital period are the same, i.e., they are
  tidally-locked, but different from the primary's spin period. In this case, only
  an orbital period is given. The lightcurve may be flat or bowed between events.
  If flat, the presumption is that the satellite is nearly spheroidal and the
  rotation is still tidally-locked to the orbit. If the lightcurve shows an overall
  "bowed" shaped, this is presumed to indicate a significantly elongated satellite.

  Fully-synchronous    (Binary Type Flag: F)  Example: Pluto/Charon, 90 Antiope
  The rotation period of the primary and satellite are the same and is the
  same as the orbital period of the satellite. In this case, the primary
  rotation period and lightcurve amplitude is given and matches the orbital
  period of the two bodies. No secondary period is given.

  Undetermined  (Binary Type Flag: U)
  Usually reserved for binaries discovered by imaging with Hubbler or very-large 
  ground-based telescopes. In most cases, the orbital parameters are not or 
  very poorly known and there are no lightcurves to determine the actual type 
  of binary, e.g., if the satellite is tidally locked to its orbital period.

In some asynchronous systems, it is not always possible to determine with
certainty which of the two is the primary and which is the secondary in such
systems. In these cases, we are forced to give the period and amplitude of one
body as that of the "primary" and the other period and amplitude as that of
the "secondary" when, in fact, the roles may be reversed from our selection.

Only the detail line of the binary report includes the type flag. It may be
followed by '?' to indicate that the object is suspected to be binary but
definitive proof was lacking.

For multiple systems and in _most_ cases, the satellite information is for the
first one discovered. In some cases, e.g., 3749 Balam, the first discovery was
for a satllite with a long orbital period of 1920 hours. It is assumed that the
satellite's rotation is not equal to the orbital period. A second satellite was
found that has a rotation period that is tidally-locked to its orbital period of
about 33.4 hours.

The Details lines also give the primary rotation period and amplitude and
secondary/orbital periods/amplitudes as appropriate. If available the estimated
effective diameter ratio (Ds/Dp) is given, as are the ratio of the semi-major
axis of the satellite orbit to the diameter of the primary (A/Dp).

The Ds/Dp ratio is a minimum in most cases since total eclipses were not seen
in the satellite's lightcurve. As of this time, there is no flag to indicate
which of the values are minimums and those that are more definitive.


4.3.1  LC_BINARY_PUB COLUMN MAPPING
-----------------------------------

Field       Header     Format      Pos   Notes
--------------------------------------------------------------------------------
Number      Number     I7          1-7   Summary: MPC assigned number. Blank if
                                         none.
EntryFlag              A1            9   * = New entry since last release
Name        Asteroid   A30       11-40   Summary: Name of asteroid
ShortRef    Reference  A30       11-40   Details: Publication reference
Type        BT         A1           42   Summary
                                           ? = Suspected
                                           B = Confirmed binary
                                           M = Multiple binary
                                         Details
					   A = Fully-asynchronous
                                           S = Singly-asynchronous
                                           F = Fully-synchronous
                                           U = Undetermined/uncertain
Status                 A1           44   
Period      Period     F13.8     46-58   Summary: Adopted period (usually
                                         primary of system)
                                         Details: Primary period
Amp1        Amp1       F4.2      60-63   Summary: Adopted maximum amplitude
                                         Details: Amplitude for primary period
Period2     Period2    F13.8     65-77   Details: Secondary period
Amp2        Amp2       F4.2      79-82   Details: Secondary amplitude
POrb        POrb       F13.8     84-96   Details: Orbital period of satellite
DsDp        DsDp       F4.2     98-101   Details: Diameter ratio primary/satellite
ADp         ADp        F6.3    103-108   Details: Ratio Oribit SMA/Primar diameter



4.4.0  LC_SPINAXIS_PUB.TXT  (POLE SOLUTIONS)
--------------------------------------------

This file includes a summary line and associated details line(s) for any
asteroid for which spin axis information has been reported. See section
2.1.1, "Spin Axis Catalogs" for additional resources and references.

The Summary line gives the number (if any) and name of the object along with the
adopted period and amplitude. The Details line gives the additional information
such as publication reference, quality rating, period, and pole solutions.


4.4.1  Q (QUALITY) RATING
-------------------------

The Q value gives our assessment of the quality of the pole solution and 
closely parallels the U rating assigned to lightcurves. 


0   Later proven to be wrong.

1   May be completely wrong. 

2  Good determination, pole likely correct to ±15-20°, but may be ambiguous with
   two or more solutions that are possible, or the sense of rotation is not 
   determined. 

3   Reliable determination of both spin axis direction and sense of rotation,
    i.e., prograde or retrograde.

4   Reliable determination. Same as Q = 3 but a shape model has been
    determined by inversion, AO, or radar, or ambiguous models have been
    eliminated by occultation profile matching.

P   A prograde rotation has been determined but no specific pole position 
    has been determined. 

R   A retrograde rotation has been determined by no specific pole position
    has been determined.

If the Q value is blank, the given pole solution has not yet been reviewed under
the new rating system.


4.4.2  LC_SPINAXIS_PUB COLUMN MAPPING
-------------------------------------

Field       Header     Format      Pos   Notes
--------------------------------------------------------------------------------
Number      Number     I7          1-7   Summary: MPC assigned number, blank if
                                         none.
EntryFlag              A1            9   * = New entry since last release.
Name        Asteroid   A30       11-40   Summary: Name or designation.
ShortRef    Reference  A30       11-40   Details: Publication reference.
Quality     Q          A3        42-44   Summary: U rating for given synodic period (see notes below)
                                         Details: Quality of pole solution (see notes below)
Period      Period     F13.8     46-58   Summary: Adopted period.
                                         Details: Period of pole solution,
                                                  usually sidereal.
AmpMax      Amp        F4.2      60-63   Summary: Adopted maximum amplitude.
                                         Details: Amplitude for given period.
L1          L1         F5.1      65-69   Ecliptic longitude of first solution.
B1          B1         F5.1      71-75   Ecliptic latitude of first solution.
L2          L2         F5.1      77-81   Ecliptic longitude of second solution.
B2          B2         F5.1      83-87   Ecliptic latitude of second solution.
L3          L3         F5.1      89-93   Ecliptic longitude of third solution.
B3          B3         F5.1      95-99   Ecliptic latitude of third solution.
L4          L4         F5.1    101-105   Ecliptic longitude of fourth solution.
B4          B4         F5.1    107-111   Ecliptic latitude of fourth solution.
SidPeriod   Sid Per    F13.8   113-125   Sidereal period of spin axis solution.
ShapeModel  Mod        A1          128   Y = Shape model reported.

NOTES
-----
The quality rating on the summary line follows the usual LCDB rating, e.g.,
2, 2+, 3-, etc. The rating on the details line follows the ratings used on the
DAMIT site, which add 0.5 for a half-step, e.g., 1.5. This was retained since
it is not possible to know if 1.5 would be 1+ or 2-.

Sometimes an entry will have a value of L1 = -1. This indicates no longitude 
was reported. The value of B1 has two interpretations

  If the latitude is -99.9, then no latitude was reported. This entry _must_ 
  have a Q value of P or R, meaning prograde or retrograde rotation was 
  determined. This is usually by seeing how the synodic rotation period changed 
  before, at, and after opposition. Other techniques than lightcurve inversion 
  can also produce a sense of rotation but no longitude/latitude pair.

  If the latitude is |b| <= 90.0, then a latitude only solution was found, 
  although it is usually more a "best guess" and can have substatial errors. 
  Again, the Q rating must be P or R, which is assigned on the basis that 
  positive latitudes imply prograde rotation and negative latitudes imply 
  retrograde rotation.  

A negative longitude will not appear for Long2-Long4.


4.5.0  LC_NPA_PUB.TXT  (NON-PRINCIPAL AXIS ROTATION - TUMBLING)
---------------------------------------------------------------

This file includes a Summary line and associated Detail line(s) for any asteroid
for which NPA (tumbling) data has been reported.

Note that the field spacing is not always 1-column wide.


4.5.1 PAR RATING
----------------

The PAR rating is adopted from Pravec et al. (2005), Icarus 173, 108-131.

Those asteroids where we believe the claims of NPA are not justified are not
included in the NPA table, nor are they given an entry in the NOTES field of
the Summary and/or Details tables. Instead a special entry is made in the
"free-form" NOTESEX field of the appropriate table and so appear in the
LC_NOTESEX file.

Following is a brief description of the PAR codes. See the Pravec paper for a
more detailed explanation.

  -4  Physical model of the NPA rotation constructed

  -3  NPA rotation reliably detected with the two periods resolved. There may
      be some ambiguities in one or both periods.

  -2  NPA rotation detected based on deviations from a single period but the
      second period is not resolved.

  -1  NPA rotation possible, i.e., deviations from a single period are seen,
      but not conclusively.

   0  Insufficient data to determine if rotation is PA or NPA

  +1  PA rotation is consistent with the data but coverage is insufficient.

  +2  PA rotation likely, or deviations from PA are small with some overlapping
      data fitting a PA rotation period.

  +3  PA rotation quite likely

  +4  PA spin vector obtained.

Entries with a positive number are rare and used when the asteroid was thought
to be tumbling but further examination showed it was likely in PA rotation, or
when the damping time to PA rotation is sufficiently long that the given
asteroid would more likely be in NPA than PA rotation.



4.5.2  LC_NPA_PUB COLUMN MAPPING
--------------------------------

Field       Header      Format    Pos   Notes
--------------------------------------------------------------------------------
Number      Number      I7        1-7   MPC assigned number. Blank if none.
EntryFlag               A1          9   * = New entry since last release.
Name        Asteroid    A30     11-40   Summary: Name or designation.
ShortRef    Reference   A30     11-40   Details: Publication reference.
Notes       PAR         A2      42-43   Summary: T flag from NOTES field*
PAR         PAR         A2      45-46   Details: NPA Probability rating.
Period      Period      F13.8   48-60   Summary: Adopted period.
                                        Details: Primary period.
AmpMax      Amp1        F4.2    62-65   Summary: Adopted maximum amplitude.
                                        Details: Amplitude for period.
Period2     Per2        F13.8   67-79   Details: Secondary period.
Amp2        Amp2        F4.2    81-84   Details: Secondary period amplitude.
                                                 

*See the discussion under the NOTES field for the Summary table regarding the
 T flag.

4.6.0  LC_REF.TXT  (REFERENCES)
-------------------------------

The LCREF.TXT file contains the complete literature citation for each reference
in the data file.  The first column gives the "shorthand" form as it appears in
the table. The second column gives the 19-character BibCode. The third column 
gives the full citation, listing all authors (up to the first four), year, 
journal, volume and page numbers.

An HTML version of the LC_REF file is also included. When loaded into a web 
browser, it becomes possible to click on the BibCode link. This takes the 
browers to the ADS web site entry for the BibCode. In many cases, it will then
be possible to download a copy of the original reference.


4.6.1  LC_REF COLUMN MAPPING
----------------------------

Field       Format     Pos     Notes
--------------------------------------------------------------------------------
EntryFlag   A1           1     * New record since last release
ShortRef    A30       3-32     Primary author and year
BibCode     A19      34-52     19-character BibCode
Citation    A1024    54-80     Full citation. This can span multiple lines. The
                               first line always starts in column 54. The second
                               and subsequent lines start in column 56 (two
                               space indent) to make the file more readable. All
                               lines go to a maximum of column 80.

N.B. Initials for names are packed, e.g., Warner, B.D. and not Warner, B. D.  


4.6.2  LC_REF ARTICLE IDS
-------------------------
Many journals are moving to giving an article ID instead of page numbers. This
makes it easier to compile the journal as papers arrive and included in a given
issue. As of mid-2011, the citation entry will indicate this where possible and
the information is available by giving volume, 'A' + article number. For example,

  Astron. J. 145, A75

Legacy entries will be updated as cicumstances allow. 


4.6.3  LC_REF WEB REFERENCES
----------------------------

Many observers choose to post their results only on a web page, or post the 
results to make them available prior to formal publication (e.g., the CALL 
web site). The ShortRef entry for references to these postings ends in "web" 
(e.g., Warner 2011web).

Prior to 2011 July, all postings on a web site were given as the current 
year + 'w'. This did not allow showing changes in the postings over the years. 
For example, an observer might post results in 2009 and then again in 2011, 
possibly being completely different. All web entries starting in 2012 will 
carry the year the specific results were posted. So, it may be possible for 
the same observer to have multiple "web" entries in the details table for 
the same object, but they may have different years, even if the results are 
the same and even if they are based on the same data.

When the results in a web posting are published, for example, a posting on 
the CALL site becomes part of a paper in the Minor Planet Bulletin, that 
specific web entry is changed to reflect the publication, i.e., the same
result does not appear as both a web and publication entry.


4.7.0  LC_COLORINDEX.TXT  (COLOR INDICES)
-----------------------------------------

The LC_COLORINDEX file contains B-V, B-R, V-R, V-I, Sloan G-R, Sloan R-I 
and/or Sloan I-Z color indices. The Summary line gives the number (if any) 
and name of the object along with the adopted period and amplitude. The 
following Details lines include the publication reference and color index 
values.

Unless most others, this report contains all entries in the Color Index 
table, even if there is an empty Summary line, i.e., no period and/or 
amplitude reported - the so-called "orphan" records.


4.7.1  LC_COLORINDEX COLUMN MAPPING
-----------------------------------

Field       Header      Format    Pos   Notes
--------------------------------------------------------------------------------
Number      Number      I7        1-7   MPC assigned number. Blank if none.
EntryFlag               A1          9   * = New entry since last release. [1]
Name        Asteroid    A30     11-40   Summary: Name or designation.
ShortRef    Reference   A30     11-40   Details: Publication reference.
Period      Period      F13.8   42-54   Summary: Adopted period.
                                        Details: Blank.
AmpMax      Amp         F4.2    56-59   Summary: Adopted maximum amplitude.
                                        Details: Blank.
B-V         BV          F6.3    61-66   Details: B-V color index.
B-R         BR          F6.3    68-73   Details: B-R color index.
V-R         VR          F6.3    75-80   Details: B-V color index.
V-I         VI          F6.3    82-87   Details: B-V color index.
SG-R        SGR         F6.3    89-94   Details: Sloan G-R color index.
SR-I        SRI         F6.3    96-101  Details: Sloan R-I color index.
SI-Z        SIZ         F6.3   103-108  Details: Sloan R-I color index.


[1] Intermediate releases may have other flags, e.g., @, #, etc. These are 
    used during the vetting process and will be cleared or set to * for the
    annual vetted release.

4.8.0  LC_NOTESEX.TXT  (EXTENDED NOTES)
---------------------------------------

The LC_NOTESEX file includes additional notes attached to either a Summary or
Details line. These notes amplify or explain information contained in the
record.

Given the nature of the file, there is no header line. Also, the formatting
of the file is slightly different depending on whether the note belongs to a
Summary or Details line. See the column mapping below.

The length of the NOTES field is not limited to a fixed-number. All notes are
forced to wrap at column 80 to the beginning of the next line. The Summary line
notes start immediately under the start of the name or designation while the
Details notes are indented in a block in order to allow the publication
reference to be included on the first line.



4.8.1  LC_NOTESEX COLUMN MAPPING
--------------------------------

Field       Format    Pos   Notes
--------------------------------------------------------------------
Number      I7        1-7   Summary: MPC assigned number. Blank if none.
                            Details: Blank.
Name        A30      9-38   Summary: Name or designation.
ShortRef    A30      9-38   Details: Publication reference.
Notes                9-80   Summary: Summary record notes.
Notes               32-80   Details: Details record notes.



4.9.0  LC_DIAMETERS.TXT  (MEASURED DIAMETERS)
-----------------------------------------

The LC_DIAMETERS file contains diameters that were measured by thermal, 
radar, occultation, or other observations. The Summary line gives the 
adopted values for period, if one is reported, H, G, Albedo, and Diameter. 
The Details line gives the H, G, Albedo, and Diameter values from the 
the publication reference along with the reported errors in these values, 
if any.

Unless most others, this report contains all entries in the Details 
table that have a diameter, even if there is an empty Summary line, i.e., 
no period and/or amplitude reported - the so-called "orphan" records.



4.9.1  LC_DIAMETERS COLUMN MAPPING
-----------------------------------

Field       Header      Format    Pos   Notes
--------------------------------------------------------------------------------
Number      Number      I7        1-7   MPC assigned number. Blank if none.
EntryFlag   <NONE>      A1          9   Entry Flag (blank or *) [1]
Name        Asteroid    A30     11-40   Summary: Name or designation.
ShortRef    Reference   A30     11-40   Details: Publication reference.
Period      Period      F13.8   42-54   Summary: Adopted period.
                                        Details: Blank.
H           H           F6.3    56-61   Summary: Adopted H.
                                        Details: H used to find diameter.
HErr        HErr        F6.3    63-68   Summary: Blank.
                                        Details: Error in H (mags)
G           G           F6.3    70-75   Summary: Adopted G.
                                        Details: G used to find diameter.
GErr        GErr        F6.3    77-83   Summary: Blank.
                                        Details: Error in G.
Albedo      Albedo      F6.4    84-89   Summary: Adopted albedo (p_V)
                                        Details: Measured/derived albedo.
AlbedoErr   A.E.        F6.4    91-96   Summary: Blank.
                                        Details: Error in albedo.
Diameter    Diam        F8.3   98-105   Summary: Adopted diameter (km)
                                        Details: Measured/derived diameter
AlbedoErr   D.E.        F8.3  107-114   Summary: Blank.
                                        Details: Error in diameter.
Notes       Notes       A5    116-120   Summary: Notes [2]
                                        Details: Notes [2]

[1] Intermediate releases may have other flags, e.g., @, #, etc. These are 
    used during the vetting process and will be cleared or set to * for the
    annual vetted release.

[2] See "Notes" in section 4.1.4



5.0.0  WIDE-FIELD SURVEYS AND SPARSE DATA SETS
--------------------------------------------------------------------------------

5.1.0  WIDE-FIELD SURVEYS
-------------------------
Since about 2012, several papers have been published that made use of the 
Palomar Transit Factory survey. The combined papers have added more than 
10,000 new asteroids to the LCDB. 

The first two papers, Polishook 2012b and Chang 2014b, produced a 
manageable number of lightcurves, meaning that each one was reviewed by 
the LCDB authors and a U code assigned. More recent papers, e.g., 
Waszczak 2015 and Chang 2015, produced thousands of "reliable" lightcurves. 
For these two papers, and those likely to follow, it is not possible to 
review each lightcurve. 

As a result, the LCDB authors adopted the policy of assigning U = 2 to any 
lightcurve from the two later papers above where the original authors claimed 
to have found a reliable period. The Waszczak paper found a number of periods 
that were not considered reliable. These were assigned U = 1 in the LCDB. Where 
no period was reported, what information that was available, e.g., amplitude, 
was entered and no U code was assigned.

Over time, at least some individual entries will be reviewed and, if necessary,
the summary line will be updated.

It should be noted that the Waszczak paper reported observations on more than 
50,000 asteroids but found less than 10,000 "reliable" periods, or about a 
16% success rate. The Chang 2015 paper had about a 27% success rate. As a 
result, statistical studies should use the wide-field data with some caution 
since they have likely introduced substantial biases, e.g., against super-fast 
or super-slow rotators, tumblers, binary objects, and, probably most 
significantly, against objects with low amplitudes, i.e., A < ~0.10-0.12 mag.

For a detailed look at these issues, the reader is referred to 

  Warner, B.D., Harris, A.W. (2011) "Using sparse photometric data sets for 
    asteroid lightcurve studies." Icarus 216, 610-624.

  Harris, A.W., Pravec, P., Warner, B.D. (2012) "Looking a gift horse in the 
    mouth: Evaluation of wide-field asteroid photometric surveys." 
    Icarus 221, 226-235.

The LCDB authors anticipate having to handle the results based on similar 
large surveys in the future. These will be handled on a case-by-case basis.


5.2.0  SPARSE DATA SETS
-----------------------
Sparse data sets differ from their wide-field counterparts in that they 
are generally the result of surveys such as the Catalina Sky Survey, i.e., 
2-5 data points a night on a few nights each lunation over several years. 
A more extensive example is the Lowell Lightcurve Database (Bowell et al., 
2014) that includes hundreds of observations for some asteroid over 10-15 
years. 

These sparse data sets can be used exclusively in shape modeling and spin 
axis studies (e.g., Durech et al., 2016). If a period and/or spin axis is 
generated by the use of sparse data _only_, the SPARSEDATA field will be 
set to TRUE in the details and/or summary record.

Depending on the number of data points involved, the result may be given 
a U rating of 2 even though the NOTES flag will include 'N' for no published
lightcurve. If the sole period result, the value is included on the summary 
line, which will also have the SPARSEDATA field set along with the U rating 
field. 

In those cases where sparse data are combined with dense lightcurves, the 
sparse data field will NOT be set but, as above, the result may be given 
a U rating and 'N' in the in NOTES field and, if the only period result, 
migrated to the summary line.

