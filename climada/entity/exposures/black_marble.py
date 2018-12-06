"""
This file is part of CLIMADA.

Copyright (C) 2017 CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---

Define BlackMarble class.
"""

__all__ = ['BlackMarble']

import logging
import math
import numpy as np
from numpy.polynomial.polynomial import polyval
from scipy import ndimage
import shapely.vectorized
from cartopy.io import shapereader
from iso3166 import countries as iso_cntry

from climada.entity.exposures.base import Exposures
from climada.util.constants import SYSTEM_DIR, ONE_LAT_KM
from climada.util.coordinates import coord_on_land
from climada.entity.exposures import nightlight as nl_utils
from climada.util.finance import gdp, income_group

LOGGER = logging.getLogger(__name__)

DEF_RES_NOAA_KM = 1
""" Default approximate resolution for NOAA NGDC nightlights in km."""

DEF_RES_NASA_KM = 0.5
""" Default approximate resolution for NASA's nightlights in km."""

DEF_HAZ_TYPE = 'TC'
""" Default hazard type used in impact functions id. """

DEF_POLY_VAL = [0, 0, 1]
""" Default polynomial transformation used. """

class BlackMarble(Exposures):
    """Defines exposures from night light intensity, GDP and income group.
    Attribute region_id is defined as:
    - United Nations Statistics Division (UNSD) 3-digit equivalent numeric code
    - 0 if country not found in UNSD.
    - -1 for water
    """

    def __init__(self):
        """ Empty initializer. """
        Exposures.__init__(self)

    def set_countries(self, countries, ref_year=2016, res_km=None,
                      sea_res=(0, 1), from_hr=None, **kwargs):
        """ Model countries using values at reference year. If GDP or income
        group not available for that year, consider the value of the closest
        available year.

        Parameters:
            countries (list or dict): list of country names (admin0) or dict
                with key = admin0 name and value = [admin1 names]
            ref_year (int, optional): reference year. Default: 2016
            res_km (float, optional): approx resolution in km. Default:
                nightlights resolution.
            sea_res (tuple, optional): (sea_coast_km, sea_res_km), where first
                parameter is distance from coast to fill with water and second
                parameter is resolution between sea points.
            from_hr (bool, optional): force to use higher resolution image,
                independently of its year of acquisition.
            kwargs (optional): 'gdp' and 'inc_grp' dictionaries with keys the
                country ISO_alpha3 code. 'poly_val' polynomial transformation
                [1,x,x^2,...] to apply to nightlight (DEF_POLY_VAL used if not
                provided). If provided, these are used.
        """
        self.clear()
        shp_file = shapereader.natural_earth(resolution='10m',
                                             category='cultural',
                                             name='admin_0_countries')
        shp_file = shapereader.Reader(shp_file)

        cntry_info, cntry_admin1 = country_iso_geom(countries, shp_file)
        fill_econ_indicators(ref_year, cntry_info, shp_file, **kwargs)

        nightlight, coord_nl, fn_nl, res_fact, res_km = get_nightlight(\
            ref_year, cntry_info, res_km, from_hr)

        # TODO parallel over _set_one_country
        for cntry_iso, cntry_val in cntry_info.items():
            LOGGER.info('Processing country %s.', cntry_val[1])
            self.append(self._set_one_country(cntry_val, nightlight, \
                coord_nl, fn_nl, res_fact, res_km, cntry_admin1[cntry_iso], \
                **kwargs))

        self.add_sea(sea_res)

    def select(self, reg_id):
        """ Select exposures with input region.

        Parameters:
            reg_id (int, str): integer iso equivalent country numeric code or
                string iso alpha-3 or alpha-2 code or country name.

        Returns:
            Exposures
        """
        if isinstance(reg_id, int):
            return Exposures.select(self, reg_id)

        try:
            return Exposures.select(self, \
                int(iso_cntry.get(reg_id).numeric))
        except KeyError:
            LOGGER.info('No country %s found.', reg_id)
            return None

    def add_sea(self, sea_res):
        """ Add sea to geometry's surroundings with given resolution.

        Parameters:
            sea_res (tuple): (sea_coast_km, sea_res_km), where first parameter
                is distance from coast to fill with water and second parameter
                is resolution between sea points
        """
        if sea_res[0] == 0:
            return

        LOGGER.info("Adding sea at %s km resolution and %s km distance from coast.",
                    str(sea_res[1]), str(sea_res[0]))

        sea_res = (sea_res[0]/ONE_LAT_KM, sea_res[1]/ONE_LAT_KM)

        min_lat = max(-90, float(np.min(self.coord.lat)) - sea_res[0])
        max_lat = min(90, float(np.max(self.coord.lat)) + sea_res[0])
        min_lon = max(-180, float(np.min(self.coord.lon)) - sea_res[0])
        max_lon = min(180, float(np.max(self.coord.lon)) + sea_res[0])

        lat_arr = np.arange(min_lat, max_lat+sea_res[1], sea_res[1])
        lon_arr = np.arange(min_lon, max_lon+sea_res[1], sea_res[1])

        lon_mgrid, lat_mgrid = np.meshgrid(lon_arr, lat_arr)
        lon_mgrid, lat_mgrid = lon_mgrid.ravel(), lat_mgrid.ravel()
        on_land = np.logical_not(coord_on_land(lat_mgrid, lon_mgrid))

        self.coord = np.array([np.append(self.coord.lat, lat_mgrid[on_land]), \
            np.append(self.coord.lon, lon_mgrid[on_land])]).transpose()
        self.value = np.append(self.value, lat_mgrid[on_land]*0)
        self.id = np.arange(1, self.value.size+1)
        self.region_id = np.append(self.region_id,
                                   lat_mgrid[on_land].astype(int)*0 - 1)
        self.impact_id = {DEF_HAZ_TYPE: np.ones(self.value.size, int)}

    @staticmethod
    def _set_one_country(cntry_info, nightlight, coord_nl, fn_nl, res_fact,
                         res_km, admin1_geom, **kwargs):
        """ Model one country.

        Parameters:
            cntry_info (lsit): [cntry_id, cnytry_name, cntry_geometry,
                ref_year, gdp, income_group]
            nightlight (np.array): nightlight in 30arcsec ~ 1km resolution.
                Row latitudes, col longitudes
            coord_nl (np.array): nightlight coordinates: [[min_lat, lat_step],
                [min_lon, lon_step]]
            fn_nl (str): file name of considered nightlight with path
            res_fact (float): resampling factor
            res_km (float): wished resolution in km
            admin1_geom (list): list of admin1 geometries to filter
            poly_val (list): list of polynomial coefficients to apply to
                nightlight
        """
        if 'poly_val' in kwargs:
            poly_val = kwargs['poly_val']
        else:
            poly_val = DEF_POLY_VAL
        geom = cntry_info[2]

        nightlight_reg, lat_reg, lon_reg, on_land = _cut_country(geom, \
            nightlight, coord_nl)
        nightlight_reg = _set_econ_indicators(nightlight_reg, cntry_info[4], \
                                              cntry_info[5], poly_val)
        if admin1_geom:
            nightlight_reg, lat_reg, lon_reg, geom, on_land = _cut_admin1( \
                nightlight_reg, lat_reg, lon_reg, admin1_geom, coord_nl, on_land)

        LOGGER.info('Generating resolution of approx %s km.', res_km)
        nightlight_reg, lat_reg, lon_reg = _resample_land(geom, nightlight_reg,\
            lat_reg, lon_reg, res_fact, on_land)

        exp_bkmrb = BlackMarble()
        exp_bkmrb.value = np.asarray(nightlight_reg).reshape(-1,)
        exp_bkmrb.coord = np.empty((exp_bkmrb.value.size, 2))
        exp_bkmrb.coord[:, 0] = lat_reg
        exp_bkmrb.coord[:, 1] = lon_reg
        exp_bkmrb.id = np.arange(1, exp_bkmrb.value.size+1)
        exp_bkmrb.region_id = np.ones(exp_bkmrb.value.shape, int)*cntry_info[0]
        exp_bkmrb.impact_id = {DEF_HAZ_TYPE: np.ones(exp_bkmrb.value.size, int)}
        exp_bkmrb.ref_year = cntry_info[3]
        exp_bkmrb.tag.description = ("{} {:d} GDP: {:.3e} income group: {:d}"+\
            "\n").format(cntry_info[1], cntry_info[3], \
            cntry_info[4], cntry_info[5])
        exp_bkmrb.tag.file_name = fn_nl
        exp_bkmrb.value_unit = 'USD'

        return exp_bkmrb

def country_iso_geom(countries, shp_file):
    """ Get country ISO alpha_3, country id (defined as the United Nations
    Statistics Division (UNSD) 3-digit equivalent numeric codes and 0 if
    country not found) and country's geometry shape.

    Parameters:
        countries (list or dict): list of country names (admin0) or dict
            with key = admin0 name and value = [admin1 names]
        shp_file (cartopy.io.shapereader.Reader): shape file

    Returns:
        cntry_info (dict): key = ISO alpha_3 country, value = [country id,
            country name, country geometry],
        cntry_admin1 (dict): key = ISO alpha_3 country, value = [admin1
            geometries]
    """
    countries_shp = {}
    list_records = list(shp_file.records())
    for info_idx, info in enumerate(list_records):
        countries_shp[info.attributes['ADMIN'].title()] = info_idx

    cntry_info = dict()
    cntry_admin1 = dict()
    if isinstance(countries, list):
        countries = {cntry: [] for cntry in countries}
        admin1_rec = list()
    else:
        admin1_rec = shapereader.natural_earth(resolution='10m',
                                               category='cultural',
                                               name='admin_1_states_provinces')
        admin1_rec = shapereader.Reader(admin1_rec)
        admin1_rec = list(admin1_rec.records())

    for country_name, prov_list in countries.items():
        country_idx = countries_shp.get(country_name.title())
        if country_idx is None:
            options = [country_opt for country_opt in countries_shp
                       if country_name.title() in country_opt]
            if not options:
                options = list(countries_shp.keys())
            LOGGER.error('Country %s not found. Possible options: %s',
                         country_name, options)
            raise ValueError
        iso3 = list_records[country_idx].attributes['ADM0_A3']
        try:
            cntry_id = int(iso_cntry.get(iso3).numeric)
        except KeyError:
            cntry_id = 0
        cntry_info[iso3] = [cntry_id, country_name.title(),
                            list_records[country_idx].geometry]
        cntry_admin1[iso3] = _fill_admin1_geom(iso3, admin1_rec, prov_list)

    return cntry_info, cntry_admin1

def fill_econ_indicators(ref_year, cntry_info, shp_file, **kwargs):
    """ Get GDP and income group per country in reference year, or it closest
    one. Source: world bank. Natural earth repository used when missing data.
    Modifies country info with values [country id, country name,
    country geometry, ref_year, gdp, income_group].

    Parameters:
        ref_year (int): reference year
        cntry_info (dict): key = ISO alpha_3 country, value = [country id,
            country name, country geometry]
        kwargs (optional): 'gdp' and 'inc_grp' dictionaries with keys the
            country ISO_alpha3 code. If provided, these are used
    """
    for cntry_iso, cntry_val in cntry_info.items():
        cntry_val.append(ref_year)
        if 'gdp' in kwargs and kwargs['gdp'][cntry_iso] != '':
            gdp_val = kwargs['gdp'][cntry_iso]
        else:
            _, gdp_val = gdp(cntry_iso, ref_year, shp_file)
        cntry_val.append(gdp_val)

        if 'inc_grp' in kwargs and kwargs['inc_grp'][cntry_iso] != '':
            inc_grp = kwargs['inc_grp'][cntry_iso]
        else:
            _, inc_grp = income_group(cntry_iso, ref_year, shp_file)
        cntry_val.append(inc_grp)

def get_nightlight(ref_year, cntry_info, res_km=None, from_hr=None):
    """ Obtain nightlight from different sources depending on reference year.
    Compute resolution factor used at resampling depending on source.

    Parameters:
        ref_year (int): reference year
        cntry_info (dict): key = ISO alpha_3 country, value = [country id,
            country name, country geometry]
        res_km (float): approx resolution in km.
        from_hr (bool, optional):
    Returns:
        nightlight (sparse.csr_matrix), coord_nl (np.array), fn_nl (str),
        res_fact (float)
    """
    if from_hr is None and ref_year > 2013:
        from_hr = True
    elif from_hr is None and ref_year <= 2013:
        from_hr = False

    if from_hr:
        if not res_km:
            res_km = 0.5
        nl_year = ref_year
        if ref_year > 2013:
            nl_year = 2016
        else:
            nl_year = 2012
        LOGGER.info("Nightlights from NASA's earth observatory for year %s.",
                    str(nl_year))
        res_fact = DEF_RES_NASA_KM/res_km
        geom = [info[2] for info in cntry_info.values()]
        geom = shapely.ops.cascaded_union(geom)
        req_files = nl_utils.check_required_nl_files(geom.bounds)
        files_exist, _ = nl_utils.check_nl_local_file_exists(req_files, \
            SYSTEM_DIR, nl_year)
        nl_utils.download_nl_files(req_files, files_exist, SYSTEM_DIR, nl_year)
        # nightlight intensity with 15 arcsec resolution
        nightlight, coord_nl = nl_utils.load_nightlight_nasa(geom.bounds, \
            req_files, nl_year)
        fn_nl = [file.replace('*', str(nl_year)) for idx, file \
                 in enumerate(nl_utils.BM_FILENAMES) if req_files[idx]]
        fn_nl = ' + '.join(fn_nl)
    else:
        if not res_km:
            res_km = 1.0
        nl_year = ref_year
        if ref_year < 1992:
            nl_year = 1992
        elif ref_year > 2013:
            nl_year = 2013
        LOGGER.info("Nightlights from NOAA's earth observation group for year %s.",
                    str(nl_year))
        res_fact = DEF_RES_NOAA_KM/res_km
        # nightlight intensity with 30 arcsec resolution
        nightlight, coord_nl, fn_nl = nl_utils.load_nightlight_noaa(nl_year)

    return nightlight, coord_nl, fn_nl, res_fact, res_km

def _fill_admin1_geom(iso3, admin1_rec, prov_list):
    """Get admin1 polygons for each input province of country iso3.

    Parameters:
        iso3 (str): admin0 country name in alpha3
        admin1_rec (list): list of admin1 records
        prov_list (list): province names
    Returns:
        list(geometry)
    """
    prov_geom = list()

    for prov in prov_list:
        found = False
        for rec in admin1_rec:
            if prov == rec.attributes['name'] and \
            rec.attributes['adm0_a3'] == iso3:
                found = True
                prov_geom.append(rec.geometry)
                break
        if not found:
            options = [rec.attributes['name'] for rec in admin1_rec \
                       if rec.attributes['adm0_a3'] == iso3]
            LOGGER.error('%s not found. Possible provinces of %s are: %s',
                         prov, iso3, options)
            raise ValueError

    return prov_geom

def _cut_admin1(nightlight, lat, lon, admin1_geom, coord_nl, on_land):
    """Cut nightlight image on box containing all the admin1 territories.

    Parameters:
        nightlight (np.array): nightlight values
        lat (np.array): latitude values in meshgrid
        lon (np.array): longitude values in meshgrid
        admin1_geom (list(shapely.geometry)): all admin1 geometries
        coord_nl (np.array): nightlight coordinates: [[min_lat, lat_step],
            [min_lon, lon_step]]
        on_land (np.array): array with true values in land points. same size
            as nightlight, lat, lon

    Returns:
        nightlight_reg, lat_reg, lon_reg (2d arrays with nightlight values,
        and coordinates in a square containing the admin1)
        on_land_reg (2d array of same size as previous with True values on land
        points)


    """
    all_geom = shapely.ops.cascaded_union(admin1_geom)

    in_lat = math.floor((all_geom.bounds[1] - lat[0, 0])/coord_nl[0, 1]), \
             math.ceil((all_geom.bounds[3] - lat[0, 0])/coord_nl[0, 1])
    in_lon = math.floor((all_geom.bounds[0] - lon[0, 0])/coord_nl[1, 1]), \
             math.ceil((all_geom.bounds[2] - lon[0, 0])/coord_nl[1, 1])

    nightlight_reg = nightlight[in_lat[0]:in_lat[-1]+1, :] \
                               [:, in_lon[0]:in_lon[-1]+1]
    nightlight_reg[nightlight_reg < 0.0] = 0.0

    lat_reg, lon_reg = np.mgrid[lat[0, 0] + in_lat[0]*coord_nl[0, 1]:
                                lat[0, 0] + in_lat[1]*coord_nl[0, 1]:
                                complex(0, nightlight_reg.shape[0]),
                                lon[0, 0] + in_lon[0]*coord_nl[1, 1]:
                                lon[0, 0] + in_lon[1]*coord_nl[1, 1]:
                                complex(0, nightlight_reg.shape[1])]

    on_land_reg = on_land[in_lat[0]:in_lat[-1]+1, :] \
                         [:, in_lon[0]:in_lon[-1]+1]

    return nightlight_reg, lat_reg, lon_reg, all_geom, on_land_reg

def _cut_country(geom, nightlight, coord_nl):
    """Cut nightlight image on box containing all the land.

    Parameters:
        geom (shapely.geometry): geometry of the region to consider
        nightlight (sparse.csr_matrix): nightlight values
        coord_nl (np.array): nightlight coordinates: [[min_lat, lat_step],
            [min_lon, lon_step]]

    Returns:
        nightlight_reg, lat_reg, lon_reg (2d arrays with nightlight values,
        and coordinates in a square containing the country)
        on_land_reg (2d array of same size as previous with True values on land
        points)
    """
    in_lat = math.floor((geom.bounds[1] - coord_nl[0, 0])/coord_nl[0, 1]), \
             math.ceil((geom.bounds[3] - coord_nl[0, 0])/coord_nl[0, 1])
    in_lon = math.floor((geom.bounds[0] - coord_nl[1, 0])/coord_nl[1, 1]), \
             math.ceil((geom.bounds[2] - coord_nl[1, 0])/coord_nl[1, 1])

    nightlight_reg = nightlight[in_lat[0]:in_lat[1]+1, in_lon[0]:in_lon[-1]+1].\
        todense()
    lat_reg, lon_reg = np.mgrid[coord_nl[0, 0] + in_lat[0]*coord_nl[0, 1]:
                                coord_nl[0, 0] + in_lat[1]*coord_nl[0, 1]:
                                complex(0, nightlight_reg.shape[0]),
                                coord_nl[1, 0] + in_lon[0]*coord_nl[1, 1]:
                                coord_nl[1, 0] + in_lon[1]*coord_nl[1, 1]:
                                complex(0, nightlight_reg.shape[1])]

    on_land_reg = np.zeros(lat_reg.shape, bool)
    for poly in geom:
        in_lat = math.floor((poly.bounds[1] - lat_reg[0, 0])/coord_nl[0, 1]), \
                 math.ceil((poly.bounds[3] - lat_reg[0, 0])/coord_nl[0, 1])
        in_lon = math.floor((poly.bounds[0] - lon_reg[0, 0])/coord_nl[1, 1]), \
                 math.ceil((poly.bounds[2] - lon_reg[0, 0])/coord_nl[1, 1])
        on_land_reg[in_lat[0]:in_lat[1]+1, in_lon[0]:in_lon[1]+1] = \
            np.logical_or( \
            on_land_reg[in_lat[0]:in_lat[1]+1, in_lon[0]:in_lon[1]+1], \
            shapely.vectorized.contains(poly, \
            lon_reg[in_lat[0]:in_lat[1]+1, in_lon[0]:in_lon[1]+1], \
            lat_reg[in_lat[0]:in_lat[1]+1, in_lon[0]:in_lon[1]+1]))

    # put zero values outside country
    nightlight_reg[np.logical_not(on_land_reg)] = 0.0

    return nightlight_reg, lat_reg, lon_reg, on_land_reg

def _resample_land(geom, nightlight, lat, lon, res_fact, on_land):
    """Model land exposures from nightlight intensities and normalized
    to GDP * (income_group + 1).

    Parameters:
        geom (shapely.geometry): geometry of the region to consider
        nightlight (np.array): nightlight values
        lat (np.array): latitude values in meshgrid
        lon (np.array): longitude values in meshgrid
        res_fact (float): resampling factor
        on_land (np.array): array with true values in land points. same size
            as nightlight, lat, lon

    Returns:
        nightlight_res, lat_res, lon_res (1d arrays with nightlight on land
        values and coordinates)
    """
    nightlight_res, lat_res, lon_res = nightlight, lat, lon
    if res_fact != 1.0:
        sum_val = nightlight.sum()
        nightlight_res = ndimage.zoom(nightlight, res_fact, mode='nearest')
        nightlight_res[nightlight_res < 0.0] = 0.0

        lat_res, lon_res = np.mgrid[
            lat[0, 0] : lat[-1, 0] : complex(0, nightlight_res.shape[0]),
            lon[0, 0] : lon[0, -1] : complex(0, nightlight_res.shape[1])]

        on_land = shapely.vectorized.contains(geom, lon_res, lat_res)

        nightlight_res[np.logical_not(on_land)] = 0.0
        nightlight_res = nightlight_res/nightlight_res.sum()*sum_val

    return nightlight_res[on_land].ravel(), lat_res[on_land], lon_res[on_land]

def _set_econ_indicators(nightlight, gdp_val, inc_grp, poly_val):
    """Model land exposures from nightlight intensities and normalized
    to GDP * (income_group + 1).

    Parameters:
        nightlight (np.matrix): nightlight values
        gdp (float): GDP to interpolate in the region
        inc_grp (float): index to weight exposures in the region
        poly_val (list): list of polynomial coefficients to apply to nightlight

    Returns:
        np.array
    """
    if nightlight.sum() > 0:
        nightlight = polyval(np.asarray(nightlight), poly_val)
        nightlight = nightlight/nightlight.sum() * gdp_val * (inc_grp+1)

    return nightlight
