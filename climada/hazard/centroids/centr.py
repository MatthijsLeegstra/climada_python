"""
This file is part of CLIMADA.

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---

Define Centroids class.
"""

import copy
import logging
from pathlib import Path
from typing import Optional, Dict, Any
import warnings

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj.crs.crs import CRS
from rasterio.warp import Resampling
from shapely.geometry.point import Point

from climada.util.constants import (DEF_CRS,
                                    ONE_LAT_KM,
                                    NATEARTH_CENTROIDS)
import climada.util.coordinates as u_coord
import climada.util.hdf5_handler as u_hdf5
import climada.util.plot as u_plot

__all__ = ['Centroids']

PROJ_CEA = CRS.from_user_input({'proj': 'cea'})

DEF_VAR_MAT = {
    'field_names': ['centroids', 'hazard'],
    'var_name': {
        'lat': 'lat',
        'lon': 'lon',
        'dist_coast': 'distance2coast_km',
        'admin0_name': 'admin0_name',
        'admin0_iso3': 'admin0_ISO3',
        'comment': 'comment',
        'region_id': 'NatId'
    }
}
"""MATLAB variable names"""

DEF_VAR_EXCEL = {
    'sheet_name': 'centroids',
    'col_name': {
        'region_id': 'region_id',
        'lat': 'latitude',
        'lon': 'longitude',
        'on_land': 'on_land',
        'dist_coast': 'dist_coast',
        'elevation': 'elevation',
       ' area_pixel': 'area_pixel'
    }
}
"""Excel variable names"""

LOGGER = logging.getLogger(__name__)


class Centroids():
    """Contains raster or vector centroids.

    Attributes
    ----------
    lat : np.array
        latitudes
    lon : np.array
        longitudes
    crs : str, optional
        coordinate reference system, default is WGS84
    area_pixel : np.array, optional
        areas
    dist_coast : np.array, optional
        distances to coast
    on_land : np.array, optional
        on land (True) and on sea (False)
    region_id : np.array, optional
        region numeric codes
    elevation : np.array, optional
        elevations
    kwargs: dicts of np.arrays, optional
        any further desired properties of centroids. Is passed to the
        GeoDataFrame constructor
    """

    def __init__(
        self,
        longitude: np.ndarray,
        latitude: np.ndarray,
        crs: str = DEF_CRS,
        region_id: Optional[np.ndarray] = None,
        on_land: Optional[np.ndarray] = None,
        dist_coast: Optional[np.ndarray] = None,
        elevation: Optional[np.ndarray] = None,
        area_pixel: Optional[np.ndarray] = None,
        **kwargs
    ):
        """Initialization

        Parameters
        ----------
        lat : np.array
            latitude of size size. Defaults to empty array
        lon : np.array
            longitude of size size. Defaults to empty array
        crs : str
            coordinate reference system
        area_pixel : np.array, optional
            area of size size. Defaults to empty array
        on_land : np.array, optional
            on land (True) and on sea (False) of size size. Defaults to empty array
        region_id : np.array, optional
            country region code of size size, Defaults to empty array
        elevation : np.array, optional
            elevation of size size. Defaults to empty array
        dist_coast : np.array, optional
            distance to coast of size size. Defaults to empty array
        """
        attr_dict = {
            'geometry': gpd.points_from_xy(longitude, latitude, crs=crs),
        }
        if region_id is not None:
            attr_dict['region_id'] = region_id
        if on_land is not None:
            attr_dict['on_land'] = on_land
        if dist_coast is not None:
            attr_dict['dist_coast'] = dist_coast
        if elevation is not None:
            attr_dict['elevation'] = elevation
        if area_pixel is not None:
            attr_dict['area_pixel'] = area_pixel
        if kwargs:
            attr_dict = dict(**attr_dict, **kwargs)
        self.gdf = gpd.GeoDataFrame(data=attr_dict, crs=crs)

    @property
    def lat(self):
        return self.gdf.geometry.y.values

    @property
    def lon(self):
        return self.gdf.geometry.x.values

    @property
    def geometry(self):
        return self.gdf['geometry']

    @property
    def on_land(self):
        return self.gdf['on_land']

    @property
    def region_id(self):
        return self.gdf['region_id']

    @property
    def elevation(self):
        return self.gdf['elevation']

    @property
    def area_pixel(self):
        return self.gdf['area_pixel']

    @property
    def dist_coast(self):
        return self.gdf['dist_coast']

    @property
    def crs(self):
        return self.gdf.crs

    @property
    def size(self):
        return self.gdf.shape[0]

    @property
    def shape(self):
        """Get shape assuming rastered data."""
        return (np.unique(self.lat).size, np.unique(self.lon).size)

    @property
    def total_bounds(self):
        """Get total bounds (minx, miny, maxx, maxy)."""
        return self.gdf.total_bounds

    @property
    def coord(self):
        """Get [lat, lon] array."""
        return np.stack([self.lat, self.lon], axis=1)

    def __eq__(self, other):
        """Return True if two centroids equal, False otherwise

        Parameters
        ----------
        other : Centroids
            centroids to compare

        Returns
        -------
        eq : bool
        """
        return self.gdf.equals(other.gdf) & u_coord.equal_crs(self.crs, other.crs)

    @classmethod
    def from_geodataframe(cls, gdf):
        return cls(longitude=gdf.geometry.x.values, latitude=gdf.geometry.y.values, crs=gdf.crs, **gdf.drop(columns=['geometry']).to_dict(orient='list'))

    @classmethod
    def from_pnt_bounds(cls, points_bounds, res, crs=DEF_CRS):
        """Create Centroids object with meta attribute according to points border data.

        raster border = point border + res/2

        Parameters
        ----------
        points_bounds : tuple
            points' lon_min, lat_min, lon_max, lat_max
        res : float
            desired resolution in same units as points_bounds
        crs : dict() or rasterio.crs.CRS, optional
            CRS. Default: DEF_CRS

        Returns
        -------
        centr : Centroids
            Centroids with meta according to given points border data.
        """
        rows, cols, ras_trans = u_coord.pts_to_raster_meta(points_bounds, (res, -res))
        x_grid, y_grid = u_coord.raster_to_meshgrid(ras_trans, cols, rows)
        return cls(lat=y_grid, lon=x_grid, crs=crs)

    def append(self, centr):
        """Append centroids points.

        Parameters
        ----------
        centr : Centroids
            Centroids to append. The centroids need to have the same CRS.

        Raises
        ------
        ValueError

        See Also
        --------
        union : Union of Centroid objects.
        """
        if not u_coord.equal_crs(self.crs, centr.crs):
            raise ValueError(
                "The centroids have different Coordinate-Reference-Systems (CRS)")
        self.gdf = pd.concat([self.gdf, centr.gdf])

    def union(self, *others):
        """
        Create the union of centroids from the inputs.
        The centroids are combined together point by point.
        All centroids must have the same CRS.

        Parameters
        ----------
        others : any number of climada.hazard.Centroids()
            Centroids to form the union with

        Returns
        -------
        centroids : Centroids
            Centroids containing the union of the centroids in others.

        """
        centroids = copy.deepcopy(self)
        for cent in others:
            centroids.append(cent)

        # remove duplicate points
        return Centroids.remove_duplicate_points(centroids)

    def get_closest_point(self, x_lon, y_lat, scheduler=None):
        """Returns closest centroid and its index to a given point.

        Parameters
        ----------
        x_lon : float
            x coord (lon)
        y_lat : float
            y coord (lat)
        scheduler : str
            used for dask map_partitions. “threads”, “synchronous” or “processes”

        Returns
        -------
        x_close : float
            x-coordinate (longitude) of closest centroid.
        y_close : float
            y-coordinate (latitude) of closest centroids.
        idx_close : int
            Index of centroid in internal ordering of centroids.
        """
        close_idx = self.geometry.distance(Point(x_lon, y_lat)).values.argmin()
        return self.lon[close_idx], self.lat[close_idx], close_idx

    def set_region_id(self, scheduler=None):
        """Set region_id as country ISO numeric code attribute for every pixel or point.

        Parameters
        ----------
        scheduler : str
            used for dask map_partitions. “threads”, “synchronous” or “processes”
        """
        ne_geom = self._ne_crs_geom(scheduler)
        LOGGER.debug('Setting region_id %s points.', str(self.size))
        self.gdf.region_id = u_coord.get_country_code(
            ne_geom.geometry[:].y.values, ne_geom.geometry[:].x.values)

    def set_area_pixel(self, min_resol=1.0e-8, scheduler=None):
        """Set `area_pixel` attribute for every pixel or point (area in m*m).

        Parameters
        ----------
        min_resol : float, optional
            if centroids are points, use this minimum resolution in lat and lon. Default: 1.0e-8
        scheduler : str
            used for dask map_partitions. “threads”, “synchronous” or “processes”
        """

        res = u_coord.get_resolution(self.lat, self.lon, min_resol=min_resol)
        res = np.abs(res).min()
        LOGGER.debug('Setting area_pixel %s points.', str(self.lat.size))
        xy_pixels = self.geometry.buffer(res / 2).envelope
        if PROJ_CEA == self.geometry.crs:
            self.gdf.area_pixel = xy_pixels.area.values
        else:
            self.gdf.area_pixel = xy_pixels.to_crs(crs={'proj': 'cea'}).area.values

    def set_area_approx(self, min_resol=1.0e-8):
        """Set `area_pixel` attribute for every pixel or point (approximate area in m*m).

        Values are differentiated per latitude. Faster than `set_area_pixel`.

        Parameters
        ----------
        min_resol : float, optional
            if centroids are points, use this minimum resolution in lat and lon. Default: 1.0e-8
        """
        res_lat, res_lon = np.abs(
            u_coord.get_resolution(self.lat, self.lon, min_resol=min_resol))
        lat_unique = np.array(np.unique(self.lat))
        lon_unique_len = len(np.unique(self.lon))
        if PROJ_CEA == self.geometry.crs:
            self.gdf.area_pixel = np.repeat(res_lat * res_lon, lon_unique_len)
            return None

        LOGGER.debug('Setting area_pixel approx %s points.', str(self.size))
        res_lat = res_lat * ONE_LAT_KM * 1000
        res_lon = res_lon * ONE_LAT_KM * 1000 * np.cos(np.radians(lat_unique))
        area_approx = np.repeat(res_lat * res_lon, lon_unique_len)
        if area_approx.size == self.size:
            self.gdf.area_pixel = area_approx
        else:
            raise ValueError('Pixel area of points cannot be computed.')

    def set_elevation(self, topo_path):
        """Set elevation attribute for every pixel or point in meters.

        Parameters
        ----------
        topo_path : str
            Path to a raster file containing gridded elevation data.
        """
        self.gdf.elevation = u_coord.read_raster_sample(topo_path, self.lat, self.lon)

    def set_dist_coast(self, signed=False, precomputed=False, scheduler=None):
        """Set dist_coast attribute for every pixel or point in meters.

        Parameters
        ----------
        signed : bool
            If True, use signed distances (positive off shore and negative on land). Default: False.
        precomputed : bool
            If True, use precomputed distances (from NASA). Default: False.
        scheduler : str
            Used for dask map_partitions. "threads", "synchronous" or "processes"
        """
        if precomputed:
            self.gdf.dist_coast = u_coord.dist_to_coast_nasa(
                self.lat, self.lon, highres=True, signed=signed)
        else:
            ne_geom = self._ne_crs_geom(scheduler)
            LOGGER.debug('Computing distance to coast for %s centroids.', str(self.size))
            self.gdf.dist_coast = u_coord.dist_to_coast(ne_geom, signed=signed)

    def set_on_land(self, scheduler=None):
        """Set on_land attribute for every pixel or point.

        Parameters
        ----------
        scheduler : str
            used for dask map_partitions. “threads”, “synchronous” or “processes”
        """
        ne_geom = self._ne_crs_geom(scheduler)
        LOGGER.debug('Setting on_land %s points.', str(self.lat.size))
        self.gdf.on_land = u_coord.coord_on_land(
            ne_geom.geometry[:].y.values, ne_geom.geometry[:].x.values)

    @classmethod
    def remove_duplicate_points(cls, centroids):
        """Return a copy of centroids with removed duplicated points

        Returns
        -------
         : Centroids
            Sub-selection of centroids withtout duplicates
        """
        return cls.from_geodataframe(centroids.gdf.drop_duplicates())

    #TODO replace with nice Geodataframe util plot method.
    def plot(self, ax=None, figsize=(9, 13), **kwargs):
        """Plot centroids scatter points over earth.

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot, optional
            axis to use
        figsize: (float, float), optional
            figure size for plt.subplots
            The default is (9, 13)
        kwargs : optional
            arguments for scatter matplotlib function

        Returns
        -------
        axis : matplotlib.axes._subplots.AxesSubplot
        """
        proj_data, _ = u_plot.get_transformation(self.crs)
        proj_plot = proj_data
        if isinstance(proj_data, ccrs.PlateCarree):
            # use different projections for plot and data to shift the central lon in the plot
            xmin, ymin, xmax, ymax = u_coord.latlon_bounds(self.lat, self.lon, buffer=pad)
            proj_plot = ccrs.PlateCarree(central_longitude=0.5 * (xmin + xmax))

        if ax is None:
            ax = self.gdf.copy().to_crs(proj_plot).plot(figsize=figsize, **kwargs)
        else:
            self.gdf.copy().to_crs(proj_plot).plot(figsize=figsize, **kwargs)

        u_plot.add_shapes(ax)
        plt.tight_layout()
        return ax


    '''
    I/O methods
    '''


    @classmethod
    def from_raster_file(cls, file_name, src_crs=None, window=None,
                         geometry=None, dst_crs=None, transform=None, width=None,
                         height=None, resampling=Resampling.nearest):
        """Create a new Centroids object from a raster file

        Select region using window or geometry. Reproject input by providing
        dst_crs and/or (transform, width, height).

        Parameters
        ----------
        file_name : str
            path of the file
        src_crs : crs, optional
            source CRS. Provide it if error without it.
        window : rasterio.windows.Window, optional
            window to read
        geometry : list of shapely.geometry, optional
            consider pixels only within these shapes
        dst_crs : crs, optional
            reproject to given crs
        transform : rasterio.Affine
            affine transformation to apply
        wdith : float
            number of lons for transform
        height : float
            number of lats for transform
        resampling : rasterio.warp,.Resampling optional
            resampling function used for reprojection to dst_crs

        Returns
        -------
        centr : Centroids
            Centroids with meta attribute according to the given raster file
        """
        meta, _ = u_coord.read_raster(
            file_name, [1], src_crs, window, geometry, dst_crs,
            transform, width, height, resampling)
        lat, lon = _meta_to_lat_lon(meta)
        return cls(lon=lon, lat=lat, crs=dst_crs)


    @classmethod
    def from_vector_file(cls, file_name, dst_crs=None):
        """Create Centroids object from vector file (any format supported by fiona).

        Parameters
        ----------
        file_name : str
            vector file with format supported by fiona and 'geometry' field.
        dst_crs : crs, optional
            reproject to given crs

        Returns
        -------
        centr : Centroids
            Centroids with points according to the given vector file
        """

        centroids = cls.from_geodataframe(gpd.read_file(file_name))
        if dst_crs is not None:
            centroids.to_crs(dst_crs, inplace=True)
        return centroids


    @classmethod
    def from_excel(cls, file_name, crs):
        """Generate a new centroids object from an excel file with column names in var_names.

        Parameters
        ----------
        file_name : str
            absolute or relative file name

        Raises
        ------
        KeyError

        Returns
        -------
        centr : Centroids
            Centroids with data from the given file
        """
        df = pd.read_excel(file_name)
        centroids = cls(**df.to_dict(orient='list'), crs=crs)
        return centroids

    def write_hdf5(self, file_name):
        """Write data frame and metadata in hdf5 format

        Parameters
        ----------
        file_name : str
            (path and) file name to write to.
        """
        LOGGER.info('Writing %s', file_name)
        store = pd.HDFStore(file_name, mode='w')
        pandas_df = pd.DataFrame(self.gdf)
        for col in pandas_df.columns:
            if str(pandas_df[col].dtype) == "geometry":
                pandas_df[col] = np.asarray(self.gdf[col])

        # Avoid pandas PerformanceWarning when writing HDF5 data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)
            # Write dataframe
            store.put('centroids', pandas_df)

        store.get_storer('centroids').attrs.metadata = {'crs': CRS.from_user_input(self.crs).to_wkt()}

        store.close()


    @classmethod
    def from_hdf5(cls, file_name):
        """Create a centroids object from a HDF5 file.

        Parameters
        ----------
        file_data : str or h5
            If string, path to read data. If h5 object, the datasets will be read from there.

        Returns
        -------
        centr : Centroids
            Centroids with data from the given file
        """
        if not Path(file_name).is_file():
            raise FileNotFoundError(str(file_name))

        with pd.HDFStore(file_name, mode='r') as store:
            metadata = store.get_storer('centroids').attrs.metadata
            # in previous versions of CLIMADA and/or geopandas, the CRS was stored in '_crs'/'crs'
            crs = metadata.get('crs')
            gdf = gpd.GeoDataFrame(store['centroids'], crs=crs)

        return cls.from_geodataframe(gdf)



    def _ne_crs_geom(self, scheduler=None):
        """Return `geometry` attribute in the CRS of Natural Earth.

        Parameters
        ----------
        scheduler : str
            used for dask map_partitions. “threads”, “synchronous” or “processes”

        Returns
        -------
        geo : gpd.GeoSeries
        """
        if u_coord.equal_crs(self.gdfgeometry.crs, u_coord.NE_CRS):
            return self.gdf.geometry
        return self.gdf.geometry.to_crs(u_coord.NE_CRS)

def _meta_to_lat_lon(meta):
    """Compute lat and lon of every pixel center from meta raster."""
    xgrid, ygrid = u_coord.raster_to_meshgrid(
        meta['transform'], meta['width'], meta['height']
        )
    return ygrid.flatten(), xgrid.flatten()
