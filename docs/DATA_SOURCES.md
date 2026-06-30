# Data Sources & Acknowledgments

Third-party datasets used in this repository — both the **input** meteorology
that drives the product and the **validation** data — with the citations and
acknowledgments their data-use policies require. (The MiCASA model product
itself is cited via its upstream DOI in `PROVENANCE.txt` and the per-file netCDF
attributes; the same files also record the exact ERA5 meteo paths per output.)

## ERA5 meteorology (model input)

The diurnalization step (`diurnalize-ERA5.r`) drives the sub-daily shape of GPP
and respiration with hourly **ECMWF ERA5** surface reanalysis — downward solar
radiation (`ssrd`), 2 m air temperature (`t2m`), and (for the opt-in soil-temp
driver) level-1 soil temperature (`stl1`) and soil water (`swvl1`). Unlike the
AmeriFlux data below, ERA5 **is an input to the product**, not just validation.

- **Source:** ECMWF ERA5, from the Copernicus Climate Change Service (C3S),
  preprocessed/regridded to the TM5 1° meteo grid (the `ea` surface stream under
  `$CARBONTRACKER/METEO/tm5-nc/ec/ea/...`; the FastTrack `ea_0005` tree is the
  near-real-time fallback). The per-day tree that fed each output is recorded in
  its `meteo_source_*` netCDF attributes.
- **Coverage used:** hourly, 1989–present (the 1999–2000 spin-up plus the full
  2001-onward record).

**Citations:**

- Hersbach, H., et al. (2020), *The ERA5 global reanalysis*, Quarterly Journal of
  the Royal Meteorological Society 146:1999–2049, https://doi.org/10.1002/qj.3803.
- Hersbach, H., et al. (2023), *ERA5 hourly data on single levels from 1940 to
  present*, Copernicus Climate Change Service (C3S) Climate Data Store (CDS),
  https://doi.org/10.24381/cds.adbb2d47.

**Acknowledgment / disclaimer** (per the Copernicus licence):

> Generated using Copernicus Climate Change Service information [year(s) of the
> data used]. Neither the European Commission nor ECMWF is responsible for any
> use that may be made of the Copernicus information or data it contains.

## AmeriFlux eddy-covariance towers

The respiration-driver evaluation (`fitter_diagnostics/ec_resp_driver_validation.py`
and `ec_diurnal_shape_overlay.py`; the "eddy-covariance gate" in
[V1_TO_V2_JUSTIFICATION.md](V1_TO_V2_JUSTIFICATION.md) §2 and
[DIURNALIZATION_ALTERNATIVES.md](DIURNALIZATION_ALTERNATIVES.md) §5.4) uses
**AmeriFlux BASE** half-hourly flux data from the 14 U.S. sites below.

- **Product:** AmeriFlux BASE-BADM (half-hourly `BASE_HH`), per-site versions as listed.
- **Obtained:** 2020-05-19 (download manifest date), as raw — non-gap-filled — sensors only.
- **Use:** validation diagnostic only (nighttime NEE ≈ respiration vs air/soil
  temperature); not an input to the MiCASA flux product.

**Acknowledgment** (per AmeriFlux Data Use Policy):

> Funding for the AmeriFlux data portal was provided by the U.S. Department of
> Energy Office of Science. We thank the site principal investigators and their
> teams for collecting and sharing the eddy-covariance data cited below.

**Per-site data citations** (author(s), period, site, version, DOI):

1. Sebastien Biraud, Marc Fischer, Stephen Chan, Margaret Torn (2002-) AmeriFlux BASE US-ARM ARM Southern Great Plains site - Lamont, Ver. 9-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246027
2. Andrew Richardson (2004-) AmeriFlux BASE US-Bar Bartlett Experimental Forest, Ver. 5-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246030
3. Bill Massman (2004-) AmeriFlux BASE US-GLE GLEES, Ver. 7-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246056
4. David Hollinger (1996-) AmeriFlux BASE US-Ho1 Howland Forest (main tower), Ver. 6-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246061
5. Nathaniel Brunsell (2007-) AmeriFlux BASE US-KFS Kansas Field Station, Ver. 6-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246132
6. Nathaniel Brunsell (2012-) AmeriFlux BASE US-KLS Kansas Land Institute, Ver. 1-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1498745
7. Ankur Desai (2001-) AmeriFlux BASE US-Los Lost Creek, Ver. 14-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246071
8. Jeffrey Wood, Lianhong Gu (2004-) AmeriFlux BASE US-MOz Missouri Ozark Site, Ver. 8-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246081
9. Asko Noormets (2005-2013) AmeriFlux BASE US-NC1 NC_Clearcut, Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246082
10. Asko Noormets (2005-) AmeriFlux BASE US-NC2 NC_Loblolly Plantation, Ver. 6-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246083
11. Asko Noormets (2013-) AmeriFlux BASE US-NC3 NC_Clearcut#3, Ver. 2-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1419506
12. Asko Noormets (2009-) AmeriFlux BASE US-NC4 NC_AlligatorRiver, Ver. 2-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1480314
13. Peter D. Blanken, Russel K. Monson, Sean P. Burns, David R. Bowling, Andrew A. Turnipseed (1998-) AmeriFlux BASE US-NR1 Niwot Ridge Forest (LTER NWT1), Ver. 15-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1246088
14. Dennis Baldocchi, Siyan Ma (2001-) AmeriFlux BASE US-Ton Tonzi Ranch, Ver. 11-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1245971

The author, period, and site-name strings above are taken verbatim from each
site's `DOI_CITATION` (BADM `GRP_DOI` group in `AMF_<site>_BIF_20200430.xlsx`);
the version is the `BASE_HH` file version. Re-pull the latest citation/version
from the [AmeriFlux site pages](https://ameriflux.lbl.gov/sites/site-list-and-pages/)
before any publication, in case a site has issued a newer data version.

> **Note — data-use policy vintage.** These data were downloaded 2020-05-19,
> before AmeriFlux's January 2021 move to CC-BY-4.0. If this analysis is
> published, check whether any individual site was under the legacy
> "AmeriFlux Data Use Policy" at download (some sites then requested an offer of
> co-authorship rather than citation alone).

## Methodology references

The respiration-partitioning and temperature-response methods (Reichstein et
al. 2005; Lasslop et al. 2010; Lloyd & Taylor 1994; Falge et al. 2001; Jung et
al. 2020; etc.) are cited in the reference lists of
[DIURNALIZATION_ALTERNATIVES.md](DIURNALIZATION_ALTERNATIVES.md) and
[V1_TO_V2_JUSTIFICATION.md](V1_TO_V2_JUSTIFICATION.md).
