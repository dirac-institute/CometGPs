import numpy as np

def subsample(time, flux, flux_err=None, npoints=100, kind="random",
              time_unit="days", night_length=8.0, seed=42,
              sigma=0.0, mean_flux=0.0, dnights=[3,6]):
    """
    Subsample an existing high-resolution time series (e.g. simulated)
    using different sampling patterns.

    Parameters
    ----------
    time : numpy.ndarray
        The array of time stamps

    flux : numpy.ndarray
        The array of flux values corresponding to the time
        stamps in `time`

    npoints : int
        The number of data points to sub-sample

    kind : str, one of {random | telescope}
        The type of subsampling to perform.
        There are different sampling patterns currently implemented:
            "random" : randomly sample `npoints` data points from
                       `time` and `flux1`
            "telescope": take day/night cycle into account, but randomly
                         subsample within a "night"
            "sparse" : sparse asteroid sampling: 3 x 1-minute images, separated by 21-minute intervals,
                                                 whole night, for K nights in a row and N nights per semester
            "ztf" :

    time_unit : str
        The unit of time that the `time` array is in.

    night_length : float
        The length of a night in hours for including day/night cycles

    seed : int
        Seed for the random number generator. Useful for reproducibility.

    dnights : [min_nights, max_nights]
        The minimum and maximum number of nights between consecutive observations. Used for
        the ZTF-like sampling

    Returns
    -------
    tsmall, fsmall, (ferrsmall) : the subsampled time series
    """

    if seed is not None:
        np.random.seed(seed)

    if kind == "random":
        # generate a random list of indices for sub-sampling
        idx = np.random.choice(np.arange(0, time.shape[0], 1, dtype=int),size=npoints, replace=False)

        tsmall = time[np.sort(idx)]
        fsmall = flux[np.sort(idx)]
        if flux_err is not None:
            ferrsmall = flux_err[np.sort(idx)]
        else:
            ferrsmall = np.zeros_like(fsmall) + sigma

    elif kind == "telescope":
        # if time unit is in days, convert length of night and day
        # from hours into days
        if time_unit == "days":
            night_length /= 24.0
            tseg = time[-1] - time[0] # total length of the time series
        else:
            tseg = (time[-1] - time[0]) / 24.0

        nightly_points = int(npoints/tseg)

        tstart = time[0]
        tend = tstart + night_length

        tsmall = []
        fsmall = []

        ferrsmall = []
        while tend <= time[-1]:
            min_ind = time.searchsorted(tstart)
            max_ind = time.searchsorted(tend)

            tshort = time[min_ind:max_ind]
            fshort = flux[min_ind:max_ind]


            idx = np.random.choice(np.arange(0, tshort.shape[0], 1, dtype=int),
                           size=nightly_points, replace=False)

            tsample = tshort[np.sort(idx)]
            fsample = fshort[np.sort(idx)]

            tsmall.append(tsample)
            fsmall.append(fsample)

            if flux_err is not None:
                ferrshort = flux_err[min_ind:max_ind]
                ferrsample = ferrshort[np.sort(idx)]
                ferrsmall.append(ferrsample)
            else:
                ferrsmall.append(np.zeros_like(fsample) + sigma)

            if time_unit == "days":
                tstart += 1.0
                tend += 1.0
            elif time_unit == "hours":
                tstart += 24.0
                tend += 24.0

    elif kind == "sparse":
        if time_unit == "days":
            night_length /= 24.0
            tseg = time[-1] - time[0] # total length of the time series
        else:
            tseg = (time[-1] - time[0]) / 24.0

        tstart = time[0]
        tend = tstart + night_length

        tsmall = []
        fsmall = []

        ferrsmall = []

        start_id = np.random.randint(0, 8, size=1)

        # one minute in fractions of days
        dt = 1./(60*24.)

        # start within the 24-minute interval of observing different asteroids
        tstart += start_id*(3*dt)

        while tend <= time[-1]:

            while tstart <= tend:
                for j in range(3):
                    tind = time.searchsorted(tstart)
                    tsmall.append(time[tind])
                    fsmall.append(flux[tind])

                    if flux_err is not None:
                        ferrsmall.append(flux_err[tind])
                    else:
                        ferrsmall.append(sigma)

                    tstart += dt

                tstart += 21.*dt

            if time_unit == "days":
                tstart += (1.0-night_length)
                tend += 1.0
            elif time_unit == "hours":
                tstart += (24.0 - night_length)
                tend += 24.0
    elif kind == "ztf":

        tstart = time[0]
        total_nights = time[-1] - time[0]

        tsmall = []
        fsmall = []

        nnights = int(npoints/2.)
        if nnights > total_nights:
            raise ValueError("Number of selected nights is longer than total light curve.")

        dnights = np.random.randint(dnights[0], dnights[1], size=npoints)
        dnights[1::2] = 0
        nights = np.zeros(npoints, float) + dnights
        nights = nights.cumsum()

        if np.max(nights) > time[-1]:
            print("Warning! The requested observations are spread over more nights" +
                  " than we currently have data for! Actual number of nights will be smaller!")

        # After 1 month, given field rises 2 hours earlier
        # Let's say we start @ +/- 1 hour of midnight
        start_times = (np.random.random(npoints) * 2.0 - 1.0) * 60.0
        start_times[::2] = start_times[1::2]
        dtimes = np.random.randint(30, 51, size=npoints)
        dtimes[::2] = 0
        times = nights + (start_times  + dtimes) / 60.0 /  24.0 + (nights * -2./24./30.)

        ferrsmall = []
        times += tstart

        for t in times:
            t_ind = time.searchsorted(t)
            tsmall.append(time[t_ind])
            fsmall.append(flux[t_ind])

            if flux_err is not None:
                ferrsmall.append(flux_err[t_ind])
            else:
                ferrsmall.append(sigma)

    fsmall = np.random.normal(np.hstack(fsmall), scale=np.hstack(ferrsmall)) + mean_flux

    #print('tsmall: ' + str(type(tsmall)))

    #return tsmall, fsmall, ferrsmall

    return np.hstack(tsmall), np.hstack(fsmall), np.hstack(ferrsmall)
