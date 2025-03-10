J/A+A/673/A114      Improving the open cluster census. II.         (Hunt+, 2023)
================================================================================
Improving the open cluster census.
II. An all-sky cluster catalogue with Gaia DR3.
    Hunt E.L., Reffert S.
    <Astron. Astrophys. 673, A114 (2023)>
    =2023A&A...673A.114H        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Milky Way ; Surveys ; Clusters, open ; Positional data ; Optical
Keywords: open clusters and associations: general - methods: data analysis -
          catalogs - astrometry
Abstract:
    Data from the Gaia satellite is revolutionising our understanding of
    the Milky Way. With every new data release, there is a need to update
    the census of open clusters. We aim to conduct a blind, all-sky search
    for open clusters using 729 million sources from Gaia DR3 down to
    magnitude G~20, creating a homogeneous catalogue of clusters including
    many new objects. We use the HDBSCAN algorithm to recover clusters. We
    validate our clusters using a statistical density test and a Bayesian
    convolutional neural network for colour-magnitude diagram
    classification. We infer basic astrometric parameters, ages,
    extinctions, and distances for the clusters in the catalogue. We
    recover 7200 clusters, 2420 of which are candidate new objects and
    4780 of which crossmatch to objects in the literature, including 134
    globular clusters. A more stringent cut of our catalogue contains 4114
    highly reliable clusters, 749 of which are new. Owing to the scope of
    our methodology, we are able to tentatively suggest that many of the
    clusters we are unable to detect may not be real, including 1152
    clusters from the Milky Way Star Cluster (MWSC) catalogue that should
    have been detectable in Gaia data. Our cluster membership lists
    include many new members and often include tidal tails. Our
    catalogue's distribution traces the galactic warp, the spiral arm
    structure, and the dust distribution of the Milky Way. While much of
    the content of our catalogue contains bound open and globular
    clusters, as many as a few thousand of our clusters are more
    compatible with unbound moving groups, which we will classify in an
    upcoming work. We have conducted the largest search for open clusters
    to date, producing a single homogeneous star cluster catalogue which
    we make available with this paper.

Description:
    We present five tables.
    Clusters: the main catalogue (Table 3).
    Clusters_rejected: an addendum to Table 3 containing rejected clusters
    in the Magellanic clouds or fragments of stellar streams. Not analysed
    in the paper, but included for completeness.
    Members: members of clusters in the main catalogue.
    Members_rejected: members of clusters in the addendum to Table 3.
    Crossmatches: a table of all clusters crossmatched against with IDs
    corresponding to clusters in the main catalogue (Table B.1)

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
clusters.dat     935     7167   Main catalogue (Table 3)
members.dat      649  1291929   Member stars of clusters (corrected version)
crossma.dat      207    29945   All (non-)xmatched clusters (Table B.1)
clustrej.dat     738      621   Objects rejected from main catalogue
membrej.dat      640   490907   Member stars of rejected clusters
--------------------------------------------------------------------------------

See also:
             I/355 : Gaia DR3 Part 1. Main source (Gaia Collaboration, 2022)
    J/A+A/646/A104 : Improving the open cluster census. I. (Hunt+, 2021)

Byte-by-byte Description of file: clusters.dat
--------------------------------------------------------------------------------
   Bytes Format Units    Label      Explanations
--------------------------------------------------------------------------------
   1- 20  A20   ---      Name       Main accepted cluster name (1)
  22- 25  I4    ---      ID         ? Internal cluster ID
  27-272  A246  ---      AllNames   Comma separated list of all crossmatched
                                     cluster names
     274  A1    ---      Type       [omg] Type of object (o: open cluster,
                                     m: moving group, g: globular cluster)
 276-286  F11.8 ---      CST        [3/99] Astrometric SNR
                                     (cluster significance test) (G1)
 288-293  I6    ---      N          [10/153797] Number of member stars
 295-305  F11.8 ---      CSTt       [2.9/99] Astrometric SNR within
                                     tidal radius (G1)
 307-311  I5    ---      Nt         [10/53464] Number of stars within
                                     tidal radius
 313-324  F12.8 deg      RAdeg      Right ascension of densest point
                                     (ICRS) at Ep=2016.0
 326-337  F12.8 deg      DEdeg      Declination of densest point
                                     (ICRS) at Ep=2016.0
 339-350  F12.8 deg      GLON       Galactic longitude
 352-362  E11.4 deg      GLAT       Galactic latitude
 364-374  F11.8 deg      r50        Radius containing 50% of members within
                                     the tidal radius
 376-386  F11.8 deg      rc         Core radius (approximate estimate)
 388-398  F11.8 deg      rt         Tidal radius (approximate estimate)
 400-410  F11.8 deg      rtot       Total radius
                                     (including tidal tails, coma, etc)
 412-424  F13.8 pc       r50pc      Radius containing 50% of members in parsecs
 426-438  F13.8 pc       rcpc       Core radius in parsecs
 440-452  F13.8 pc       rtpc       Tidal radius in parsecs
 454-466  F13.8 pc       rtotpc     Total radius in parsecs
 468-480  F13.8 mas/yr   pmRA       Mean proper motion in right ascension
                                     multipled by cos(dec)
 482-492  F11.8 mas/yr s_pmRA       Standard deviation of pmRA
 494-503  F10.8 mas/yr e_pmRA       Standard error of pmRA
 505-516  F12.8 mas/yr   pmDE       Mean proper motion in declination
 518-528  F11.8 mas/yr s_pmDE       Standard deviation of pmDE
 530-539  F10.8 mas/yr e_pmDE       Standard error of pmDE
 541-552  F12.8 mas      Plx        Mean parallax
 554-564  F11.8 mas    s_Plx        Standard deviation of parallax
 566-575  F10.8 mas    e_Plx        Standard error of parallax
 577-591  F15.8 pc       dist16     16th percentile of maximum likelihood
                                     distance
 593-607  F15.8 pc       dist50     50th percentile of maximum likelihood
                                     distance
 609-624  F16.8 pc       dist84     84th percentile of maximum likelihood
                                     distance
 626-630  I5    ---      Ndist      Number of stars used for distance
                                     calculation
     632  I1    ---      globalPlx  [0/1] Flag indicating clusters for which a
                                     star-by-star parallax offset correction was
                                     not possible during distance calculation
 634-649  F16.8 pc       X          X coordinate in galactocentric galactic
                                     coordinates
 651-666  F16.8 pc       Y          Y coordinate in galactocentric galactic
                                     coordinates
 668-683  F16.8 pc       Z          Z coordinate in galactocentric galactic
                                     coordinates
 685-697  F13.8 km/s     RV         ? Mean Gaia DR3 radial velocity
 699-711  F13.8 km/s   s_RV         ? Standard deviation of radial velocity
 713-725  F13.8 km/s   e_RV         ? Standard error of radial velocity
 727-730  I4    ---    n_RV         Number of member stars with a
                                     radial velocity
 732-740  E9.4  ---      CMDCl2.5   [0/1] 2.5th percentile of CMD class
 742-750  E9.4  ---      CMDCl16    [0/1] 16th percentile of CMD class
 752-761  F10.8 ---      CMDCl50    [0/1] 50th percentile of CMD class
 763-772  F10.8 ---      CMDCl84    [0/1] 84th percentile of CMD class
 774-783  F10.8 ---      CMDCl97.5  [0/1] 97.5th percentile of CMD class
 785-787  A3    ---      CMDClHuman Human-assigned CMD class
                                     (where available) (G2)
 789-798  F10.8 [yr]     logAge16   16th percentile of logarithm of cluster age
 800-809  F10.8 [yr]     logAge50   50th percentile of logarithm of cluster age
 811-821  F11.8 [yr]     logAge84   84th percentile of logarithm of cluster age
 823-831  E9.4  mag      AV16       16th percentile of V-band cluster extinction
 833-842  F10.8 mag      AV50       50th percentile of V-band cluster extinction
 844-853  F10.8 mag      AV84       84th percentile of V-band cluster extinction
 855-863  E9.4  mag      diffAV16   16th percentile of approximate V-band
                                     differential extinction
 865-874  F10.8 mag      diffAV50   50th percentile of approximate V-band
                                     differential extinction
 876-885  F10.8 mag      diffAV84   84th percentile of approximate V-band
                                     differential extinction
 887-897  F11.8 mag      MOD16      16th percentile of photometrically
                                     estimated distance modulus
 899-909  F11.8 mag      MOD50      50th percentile of photometrically
                                     estimated distance modulus
 911-921  F11.8 mag      MOD84      84th percentile of photometrically
                                     estimated distance modulus
 923-924  I2    ---      minClSize  [-1/80] HDBSCAN parameter used to construct
                                     cluster membership list (-1 for membership
                                     lists not made with HDBSCAN)
     926  I1    ---      isMerged   [0/1] Flag indicating manually merged
                                     cluster membership list
     928  I1    ---      isGMMMemb  [0/1] Flag indicating cluster membership
                                     list constructed using additional Gaussian
                                     mixture model post-processing step
 930-931  I2    ---      NXmatches  [0/23] Number of unique crossmatches to
                                     this cluster
 933-935  A3    ---      XmatchType Type of crossmatch (G3)
--------------------------------------------------------------------------------
Note (1): HSC_1 to HSC_10 for the clusters detected in this study.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: members.dat
--------------------------------------------------------------------------------
   Bytes Format Units    Label        Explanations
--------------------------------------------------------------------------------
   1- 20  A20   ---      Name         Main accepted cluster name
  22- 25  I4    ---      ID           ? Internal cluster ID
  27- 45  I19   ---      GaiaDR3      Gaia DR3 source ID
      47  I1    ---      inrt         [0/1] Flag indicating if cluster is within
                                       tidal radius (i.e. is a good member)
  49- 58  F10.8 ---      Prob         [0/1] Membership probability
  60- 71  F12.8 deg      RAdeg        Right ascension of star
                                      (ICRS) at Ep=2016.0
  73- 82  F10.8 arcsec e_RAdeg        Standard error on RA
  84- 95  F12.8 deg      DEdeg        Declination of star (ICRS) at Ep=2016.0
  97-106  F10.8 arcsec e_DEdeg        Standard error on DE
 108-119  F12.8 deg      GLON         Galactic longitude
 121-131  E11.4 deg      GLAT         Galactic latitude
 133-143  E11.4 mas/yr   pmRA         Mean proper motion in right ascension
                                       multipled by cos(dec)
 145-154  F10.8 mas/yr e_pmRA         Standard error of pmRA
 156-166  E11.4 mas/yr   pmDE         Mean proper motion in declination
 168-177  F10.8 mas/yr e_pmDE         Standard error of pmDE
 179-189  E11.4 mas      Plx          Mean parallax
 191-200  F10.8 mas    e_Plx          Standard error of parallax
 202-212  F11.8 um-1     pscol        ? Estimated pseudocolour
 214-224  F11.8 um-1   e_pscol        ? Standard error on pseudocolor
 226-236  E11.4 ---      PlxpmRACor   Correlation between parallax and pmRA
 238-248  E11.4 ---      PlxpmDECor   Correlation between parallax and pmDE
 250-260  E11.4 ---      pmRApmDECor  Correlation between pmRA and pmdE
 262-272  E11.4 ---      PlxpscolCor  ? Correlation between parallax
                                       and pseudocolour
 274-284  E11.4 ---      pmRApscolCor ? Correlation between pmRA
                                       and pseudocolour
 286-296  E11.4 ---      pmDEpscolCor ? Correlation between pmDE
                                       and pseudocolour
 298-299  I2    ---      Solved       [31/95] Gaia DR3 flag indicating which
                                       parameters solved for
 301-312  F12.8 deg      ELAT         Ecliptic latitude
 314-323  F10.7 um-1     nueff        ? Effective wavenumber of source used
                                       in astrometric solution
 325-335  F11.8 ---      RUWE         Renormalised unit weight error
 337-346  F10.8 ---      FidelityV1   Rybizki et al. (2022MNRAS.510.2597R) V1
                                       fidelity parameter
 348-366  F19.8 e-/s     FG           Mean G-band flux
 368-383  F16.8 e-/s   e_FG           Error on mean G-band flux
 385-402  F18.8 e-/s     FBP          Mean BP-band flux
 404-420  F17.8 e-/s   e_FBP          Error on mean BP-band flux
 422-439  F18.8 e-/s     FRP          Mean RP-band flux
 441-457  F17.8 e-/s   e_FRP          Error on mean RP-band flux
 459-468  F10.7 mag      Gmag         Mean G-band magnitude
 470-479  F10.7 mag      BPmag        Mean BP-band magnitude
 481-490  F10.7 mag      RPmag        Mean RP-band magnitude
 492-502  E11.4 mag      BP-RP        BP-RP colour
 504-514  E11.4 mag      BP-G         BP-G colour
 516-526  E11.4 mag      G-RP         G-RP colour
 528-536  E9.3  km/s     RV           ? Mean Gaia DR3 radial velocity
 538-549  F12.8 km/s   e_RV           ? Standard error of radial velocity
 551-554  F4.1  ---    n_RV           ? Method used to obtain the RV
 556-561  F6.1  ---    o_RV           ? Number of transits used to compute
                                       the RV
 563-567  F5.1  ---    o_RVd          ? Number of transits that underwent
                                       deblending
 569-579  F11.7 mag      GRVSmag      ? Mean Grvs magnitude
 581-591  F11.8 mag    e_GRVSmag      ? Error on mean Grvs magnitude
 593-598  F6.1  ---    o_GRVSmag      ? Number of transits used to construct
                                       Grvs magnitude
 600-611  F12.7 km/s     Vbroad       ? Spectral line broadening parameter
 613-625  F13.8 km/s   e_Vbroad       ? Error on the spectral line broadening
 627-631  F5.1  ---    o_Vbroad       ? Number of transits used to compute
                                       vbroad
 633-645  A13   ---      VarFlag      Gaia DR3 photometric variability flag
     647  I1    ---      NSS          [0/6] Flag indicating source has
                                        additional information in the Gaia DR3
                                        non single star tables
     649  I1    ---      RVS          [0/1] Flag indicating the availability of
                                       mean RVS spectrum for this source
--------------------------------------------------------------------------------

Byte-by-byte Description of file: crossma.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label         Explanations
--------------------------------------------------------------------------------
   1-  4  I4    ---     ID            ? Internal cluster ID (blank if object has
                                       no crossmatched cluster in our work)
   6- 26  A21   ---     NameLit       Name in literature catalogue
  28- 43  A16   ---     SourceCat     Source catalogue reference
  45- 72  A28   ---     TypeSourceCat Comma separated list of all source
                                       catalogue crossmatch types for this
                                       crossmatch
  74- 83  E10.4 deg     Sep           [0/2]? Separation between cluster centres
  85- 95  F11.8 ---     SepTidal      [0/0.94]? Separation in terms of largest
                                       tidal radius
  97-107  F11.8 ---     SepTidalLit   [0/1]? Separation in terms of literature
                                       tidal radius
 109-119  F11.8 ---     SepTidalData  [0/1]? Separation in terms of tidal radius
                                       in this work
 121-130  E10.4 mas/yr  pmRASep       [0/17.41]? Separation in pmRA times
                                       cos(dec) between this work and literature
 132-142  F11.8 ---     pmRASigma     ? pmRA separation in terms of uncertainty
                                       on pmRA in this work and
                                       literature added in quadrature
 144-153  E10.4 mas/yr  pmDESep       [0/11.24]? Separation in pmDE between
                                       this work and literature
 155-164  E10.4 ---     pmDESigma     ? pmDE separation in terms of uncertainty
                                       on pmRA in this work and
                                       literature added in quadrature
 166-175  E10.4 mas     PlxSep        [0/2.8]? Separation in parallax between
                                       this work and literature
 177-186  E10.4 ---     PlxSigma      ? Parallax separation in terms of
                                       uncertainty on pmRA in this work and
                                       literature added in quadrature
 188-196  E9.3  ---     maxSigma      ? Maximum value of pmRASepSigma,
                                       pmDESepSigma, and PlxSepSigma
                                       (must be less than 2.0 for a valid
                                       astrometric match)
 198-207  E10.4 ---     meanSigma     ? Mean value of pmRASepSigma,
                                       pmDESepSigma, and PlxSepSigma
--------------------------------------------------------------------------------

Byte-by-byte Description of file: clustrej.dat
--------------------------------------------------------------------------------
   Bytes Format Units    Label      Explanations
--------------------------------------------------------------------------------
   1- 15  A15   ---      Name       Main accepted cluster name
  17- 20  I4    ---      ID         [7167/7787] Internal cluster ID
  22- 69  A48   ---      AllNames   Comma separated list of all crossmatched
                                     cluster names
      71  A1    ---      Type       [omg] Type of object (o: open cluster,
                                     m: moving group, g: globular cluster)
  73- 83  F11.8 ---      CST        [3/99] Astrometric SNR
                                     (cluster significance test) (G1)
  85- 90  I6    ---      N          [10/104862] Number of member stars
  92-102  F11.8 ---      CSTt       [3/99] Astrometric SNR within tidal radius
 104-108  I5    ---      Nt         [10/40000] Number of stars within tidal
                                     radius
 110-121  F12.8 deg      RAdeg      Right ascension of densest point
                                      (ICRS) at Ep=2016.0
 123-134  F12.8 deg      DEdeg      Declination of densest point
                                     (ICRS) at Ep=2016.0
 136-147  F12.8 deg      GLON       Galactic longitude
 149-160  F12.8 deg      GLAT       Galactic latitude
 162-172  F11.8 deg      r50        Radius containing 50% of members within
                                     the tidal radius
 174-184  F11.8 deg      rc         Core radius (approximate estimate)
 186-196  F11.8 deg      rt         Tidal radius (approximate estimate)
 198-208  F11.8 deg      rtot       Total radius
                                     (including tidal tails, coma, etc)
 210-222  F13.8 pc       r50pc      Radius containing 50% of members in parsecs
 224-236  F13.8 pc       rcpc       Core radius in parsecs
 238-251  F14.8 pc       rtpc       Tidal radius in parsecs
 253-266  F14.8 pc       rtotpc     Total radius in parsecs
 268-280  F13.8 mas/yr   pmRA       Mean proper motion in right ascension
                                     multipled by cos(dec)
 282-292  F11.8 mas/yr s_pmRA       Standard deviation of pmRA
 294-303  F10.8 mas/yr e_pmRA       Standard error of pmRA
 305-316  F12.8 mas/yr   pmDE       Mean proper motion in declination
 318-328  F11.8 mas/yr s_pmDE       Standard deviation of pmDE
 330-339  F10.8 mas/yr e_pmDE       Standard error of pmDE
 341-352  F12.8 mas      Plx        Mean parallax
 354-363  F10.8 mas    s_Plx        Standard deviation of parallax
 365-374  F10.8 mas    e_Plx        Standard error of parallax
 376-390  F15.8 pc       dist16     16th percentile of maximum likelihood
                                     distance
 392-407  F16.8 pc       dist50     50th percentile of maximum likelihood
                                     distance
 409-424  F16.8 pc       dist84     84th percentile of maximum likelihood
                                     distance
 426-430  I5    ---      Ndist      [7/46664] Number of stars used for distance
                                     calculation
     432  I1    ---      globalPlx  [0/1] Flag indicating clusters for which a
                                     star-by-star parallax offset correction was
                                     not possible during distance calculation
 434-449  F16.8 pc       X          X coordinate in galactocentric galactic
                                     coordinates
 451-466  F16.8 pc       Y          Y coordinate in galactocentric galactic
                                     coordinates
 468-483  F16.8 pc       Z          Z coordinate in galactocentric galactic
                                     coordinates
 485-497  F13.8 km/s     RV         ? Mean Gaia DR3 radial velocity
 499-511  F13.8 km/s   s_RV         ? Standard deviation of radial velocity
 513-525  F13.8 km/s   e_RV         ? Standard error of radial velocity
 527-530  I4    ---    n_RV         Number of member stars with a
                                     radial velocity
 532-540  E9.4  ---      CMDCl2.5   [0/1] 2.5th percentile of CMD class
 542-550  E9.4  ---      CMDCl16    [0/1] 16th percentile of CMD class
 552-561  F10.8 ---      CMDCl50    [0/1] 50th percentile of CMD class
 563-572  F10.8 ---      CMDCl84    [0/1] 84th percentile of CMD class
 574-583  F10.8 ---      CMDCl97.5  [0/1] 97.5th percentile of CMD class
 585-587  A3    ---      CMDClHuman Human-assigned CMD class
                                     (where available) (G2)
 589-598  F10.8 [yr]     logAge16   16th percentile of logarithm of cluster age
 600-609  F10.8 [yr]     logAge50   50th percentile of logarithm of cluster age
 611-621  F11.8 [yr]     logAge84   84th percentile of logarithm of cluster age
 623-631  E9.4  mag      AV16       16th percentile of V-band cluster extinction
 633-642  F10.8 mag      AV50       50th percentile of V-band cluster extinction
 644-653  F10.8 mag      AV84       84th percentile of V-band cluster extinction
 655-663  E9.4  mag      diffAV16   16th percentile of approximate V-band
                                     differential extinction
 665-674  F10.8 mag      diffAV50   50th percentile of approximate V-band
                                     differential extinction
 676-685  F10.8 mag      diffAV84   84th percentile of approximate V-band
                                     differential extinction
 687-697  F11.8 mag      MOD16      16th percentile of photometrically
                                     estimated distance modulus
 699-709  F11.8 mag      MOD50      50th percentile of photometrically
                                     estimated distance modulus
 711-721  F11.8 mag      MOD84      84th percentile of photometrically
                                     estimated distance modulus
 723-724  I2    ---      minClSize  [-1/80] HDBSCAN parameter used to construct
                                     cluster membership list (-1 for membership
                                     lists not made with HDBSCAN)
     726  I1    ---     isMerged    [0/1] Flag indicating manually merged
                                     cluster membership list
     728  I1    ---     isGMMMemb   [0/1] Flag indicating cluster membership
                                     list constructed using additional Gaussian
                                     mixture model post-processing step
     730  I1    ---     NXmatches   [0/2] Number of unique crossmatches to this
                                     cluster
 732-734  A3    ---     XmatchType  Type of crossmatch (G3)
     736  I1    ---     isMC        [0/1] Flag indicating cluster is believed to
                                     be in Magellanic clouds
     738  I1    ---     isStream    [0/1] Flag indicating cluster is believed to
                                     be a component of a stellar stream
--------------------------------------------------------------------------------

Byte-by-byte Description of file: membrej.dat
--------------------------------------------------------------------------------
   Bytes Format Units    Label        Explanations
--------------------------------------------------------------------------------
   1- 15  A15   ---      Name         Main accepted cluster name
  17- 20  I4    ---      ID           [7167/7787] Internal cluster ID
  22- 40  I19   ---      GaiaDR3      Gaia DR3 source ID
      42  I1    ---      inrt         [0/1] Flag indicating if cluster is within
                                       tidal radius (i.e. is a good member)
  44- 53  F10.8 ---      Prob         [0/1] Membership probability
  55- 66  F12.8 deg      RAdeg        Right ascension of star
                                       (ICRS) at Ep=2016.0
  68- 77  F10.8 arcsec e_RAdeg        Standard error on RA
  79- 90  F12.8 deg      DEdeg        Declination of star (ICRS) at Ep=2016.0
  92-101  F10.8 arcsec e_DEdeg        Standard error on DE
 103-114  F12.8 deg      GLON         Galactic longitude
 116-127  F12.8 deg      GLAT         Galactic latitude
 129-139  E11.4 mas/yr   pmRA         Mean proper motion in right ascension
                                       multipled by cos(dec)
 141-150  F10.8 mas/yr e_pmRA         Standard error of pmRA
 152-162  E11.4 mas/yr   pmDE         Mean proper motion in declination
 164-173  F10.8 mas/yr e_pmDE         Standard error of pmDE
 175-185  E11.4 mas      Plx          Mean parallax
 187-196  F10.8 mas    e_Plx          Standard error of parallax
 198-208  F11.8 um-1     pscol        ? Estimated pseudocolour
 210-220  F11.8 um-1   e_pscol        ? Standard error on pseudocolor
 222-232  E11.4 ---      PlxpmRACor   Correlation between parallax and pmRA
 234-244  E11.4 ---      PlxpmDECor   Correlation between parallax and pmDE
 246-256  E11.4 ---      pmRApmDECor  Correlation between pmRA and pmDE
 258-268  E11.4 ---      PlxpscolCor  ? Correlation between parallax and
                                       pseudocolour
 270-280  E11.4 ---      pmRApscolCor ? Correlation between pmRA and
                                       pseudocolour
 282-292  E11.4 ---      pmDEpscolCor ? Correlation between pmDE and
                                       pseudocolour
 294-295  I2    ---      Solved       [31/95] Gaia DR3 flag indicating which
                                       parameters solved for
 297-307  E11.4 deg      ELAT         Ecliptic latitude
 309-318  F10.7 um-1     nueff        ? Effective wavenumber of source used in
                                       astrometric solution
 320-330  F11.8 ---      RUWE         Renormalised unit weight error
 332-341  F10.8 ---      FidelityV1   Rybizki et al. (2022MNRAS.510.2597R)
                                       V1 fidelity parameter
 343-360  F18.8 e-/s     FG           Mean G-band flux
 362-376  F15.8 e-/s   e_FG           Error on mean G-band flux
 378-394  F17.8 e-/s     FBP          Mean BP-band flux
 396-410  F15.8 e-/s   e_FBP          Error on mean BP-band flux
 412-429  F18.8 e-/s     FRP          Mean RP-band flux
 431-446  F16.8 e-/s   e_FRP          Error on mean RP-band flux
 448-457  F10.7 mag      Gmag         Mean G-band magnitude
 459-468  F10.7 mag      BPmag        Mean BP-band magnitude
 470-479  F10.7 mag      RPmag        Mean RP-band magnitude
 481-491  E11.4 mag      BP-RP        BP-RP colour
 493-503  E11.4 mag      BP-G         BP-G colour
 505-515  E11.4 mag      G-RP         G-RP colour
 517-529  F13.8 km/s     RV           ? Mean Gaia DR3 radial velocity
 531-542  F12.8 km/s   e_RV           ? Standard error of radial velocity
 544-547  F4.1  ---    n_RV           ? Method used to obtain the RV
 549-553  F5.1  ---    o_RV           ? Number of transits used to
                                       compute the RV
 555-559  F5.1  ---    o_RVd          ? Number of transits that underwent
                                       deblending
 561-571  F11.7 mag      GRVSmag      ? Mean Grvs magnitude
 573-583  F11.8 mag    e_GRVSmag      ? Error on mean Grvs magnitude
 585-589  F5.1  ---    o_GRVSmag      ? Number of transits used to construct
                                       Grvs magnitude
 591-602  F12.7 km/s     Vbroad       ? Spectral line broadening parameter
 604-616  F13.8 km/s   e_Vbroad       ? Error on the spectral line broadening
 618-622  F5.1  ---    o_Vbroad       ? Number of transits used to compute
                                       vbroad
 624-636  A13   ---      VarFlag      Gaia DR3 photometric variability flag
     638  I1    ---      NSS          [0/6] Flag indicating source has
                                       additional information in the Gaia DR3
                                       non single star tables
     640  I1    ---      RVS          [0/1] Flag indicating the availability of
                                       mean RVS spectrum for this source
--------------------------------------------------------------------------------

Global notes:
Note (G1): set to 99 for clusters with an SNR greater than approximately 38
  (due to numerical precision limit in our pipeline/scipy).
Note (G2): Human-assigned CMD class as follows:
            FP  = false positive
            FP? = false positive ?
            TP  = true positive
            TP? = true positive ?
Note (G3): Type of crossmatch as follows:
            1:1 = one to one
            1:m = one to many
            m:1 = many to one
            m:m = many to many
            blank for new object
--------------------------------------------------------------------------------

Acknowledgements:
    Emily Hunt, ehunt(at)lsw.uni-heidelberg.de

References:
    Hunt & Reffert, Paper I   2021A&A...646A.104H, Cat. J/A+A/646/A104

History:
    16-May-2023: on-line version
    20-Jul-2023: corrected members.dattable (correct declination)

================================================================================
(End)                                        Patricia Vannier [CDS]  23-Mar-2023
