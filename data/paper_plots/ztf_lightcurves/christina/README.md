# Notes on ZTF Data

 ZTF asteroid lightcurves vary not only because of the inherent rotation (typically spanning over a couple of hours) but also
 because of a changing phase angle, the angle between the sun, asteroid, and Earth. This variation happens over several months
 but is visibly present in data captured over longer spans of time. In order to correct for this phase angle change, we 
 calculate the predicted magnitude for the source (the prediction includes distance between the earth and asteroid, distance 
 between the asteroid and the sun, and the phase angle of the asteroid) and then subtract that prediction from the measured
 magnitude — this implicitly corrects for the phase angle (as well as the distance effects).
 
 This calculation for the predicted magnitude is done with - https://github.com/dirac-institute/sso_tools/blob/master/sso_tools/ztf/objectInfo.py#L222
 (and since the predicted magnitudes are in one band only, we then also calculate mean offsets per bandpass, so that the 
 resulting “corrected” magnitudes = (observed mag - predicted mag - color offset) are approximately centered on 0)
 
 Each asteroid has a \<id\>_obs.csv file with the observations and a \<id\>_phot.png file visualizing the observations
(and for each object, there will be more coming in with the phased result from gatspy multi-band 2-term LS fitting)

In the .csv files: reading the file with pandas will give you back the columns appropriately (although I forgot to drop 
the index when I wrote it, so you may end up having to drop that again, not sure).
The new columns beyond the original alerts are: obsdate,mjd,magcorrZTF,magOO,phaseangle,heliodist,geodist,velocity,magcorrOO
 
* **obsdate** = isot formatted date 
  * (useful for astropy, or trying to match against IRSA files, but probably not totally important for you)
* **mjd** = mjd formatted date 
  * (useful for fitting but quite obvious)
* **magcorrZTF** = the observed magnitude minus the ZTF prediction for the magnitude 
  * (the ZTF prediction is a simple prediction from their system, which uses astcheck) 
  * I was using this early on, but honestly it's not much use anymore. Ignore.
* **magOO** = the predicted magnitude (at this distance and phase) from OpenOrb, with an additional correction for the mean color of the asteroid in each bandpass
  * this is more interesting for you!
  * OpenOrb predicts this magnitude based on the orbit of the object, and generates the predicted V magnitude of the object. 
  * OpenOrb assumes the distance from the Earth and Sun calculated by OOrb, and then also takes into account the phase of the object using a simple HG phase formula (note that I've used G=0.15 always here). 
  * These predicted V-band magnitudes are then corrected to the bandpass of the observation, by simply taking the mean of the differences between predicted (V) mags and the observed magnitudes in each bandpass (and then adding this offset to each predicted magnitude). (edited) 
* **phaseangle** = the phase angle of the object (angle between earth-asteroid-sun)
* **heliodist** = the heliocentric distance to the asteroid (asteroid-sun) 
* **geodist** = the geocentric distance to the asteroid (earth-asteroid)
* **velocity** = the apparent velocity of the asteroid 
  * I don't use this, but it could be relevant if the object was moving fast enough
* **magcorrOO** == the magnitude of the asteroid in the alert, minus the OpenOrb predicted magnitude
  * although then it has another correction for the bandpass (since OOrb's predicted magnitude is in V and these are observations in both r and g bands, sometimes even i)
  * this is likely where you'd like to start
  * The mean values of magcorrOO in each bandpass should be close to 0.
  
## Caveats

### Outliers

I do not do any outlier rejection before outputting the csv values. There can definitely be outliers - 
I generally reject observations which are more than 3 * rms in the photometry scatter away from the mean before 
fitting the result. Not sure if GP likes to handle this differently - but there can be problems with the photometry 
such as including bright stars instead of the actual asteroid.

### Color Correction
It's definitely a first order correction I've done. After fitting for the lightcurve, you can do a better job.

### Phase curve correction

Not all asteroids actually have G=0.15 so my phase function correction could be considered another 'first order correction', and for some of the Jovian asteroids you can already see that not having this actually correct makes a difference. This is both good and bad - bad because its' another thing that ought to be fit for (even after this first-order correction) and good because learning about G tells you more about the asteroid and if it still needs to be corrected for after this first-order correction, there's more to be measured! So you may still see effects that you would like to fit for.

For more info on phase curve correction and fitting G, look at  https://www.sciencedirect.com/science/article/pii/S001910351000151X.

A further potential wrinkle on the phase curve correction: If you wanted to fit for the phase curve, you could look at what the phase function looks like and see if it's something you could add into the GP fitting. In which case, you'd likely still want to take out the distance correction separately -- so you'd take out magnitude effects related to the heliodist and geodist and fit the remaining phase function stuff separately.

### LS Fit
Lynne thinks perhaps the only periods from LS that she believes out of this lot are for:
11351, 1173, maybe 16974, 18046, 21601, maybe 23480, maybe 2759.
