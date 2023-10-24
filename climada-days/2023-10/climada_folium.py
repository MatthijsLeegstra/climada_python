import math

import folium
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import Polygon

import climada.util.coordinates as u_coord


def estimate_zoom_level(geodf):
    """crude implementation for guessing the initial zoom level
    based on the latitudinal expansion of the geo dataframe
    TODO: find a better formula
    """
    vdist = abs(geodf.latitude.max() - geodf.latitude.min())
    return 9 - int(math.log2(vdist) / 0.8)


def center(geodf):
    cy = (geodf.latitude.max() + geodf.latitude.min()) / 2
    cx = (geodf.longitude.max() + geodf.longitude.min()) / 2
    return (cy, cx)


def voronoi_polygons(country):
    from geovoronoi import voronoi_regions_from_coords

    iso3a = u_coord.country_to_iso(country, representation="alpha3")
    country_shape = u_coord.get_land_geometry([iso3a])
    def _(geodf):
        coords = geodf.loc[:,['longitude', 'latitude']].values
        voronoi_polys, voronoi_idx = voronoi_regions_from_coords(coords, country_shape)
        remap = {vi: k for k,v in voronoi_idx.items() for vi in v}
        return pd.Series([voronoi_polys[remap.get(i, 0)  # TODO: gross! fix it!
                                           ] for i in range(geodf.shape[0])])
    return _


def raster_polygons(width):
    def _(geodf):
        return geodf.apply(lambda row: Polygon((
            (row.longitude - width/2, row.latitude - width/2),
            (row.longitude - width/2, row.latitude + width/2),
            (row.longitude + width/2, row.latitude + width/2),
            (row.longitude + width/2, row.latitude - width/2),
            (row.longitude - width/2, row.latitude - width/2),
        )), axis=1)
    return _


def make_bins(dseries, n):
    distance = dseries.max() - dseries.min()
    om = int(math.log(distance, 10))
    lower = int(dseries.min() / math.pow(10, om))
    upper = (int(dseries.max() / math.pow(10, om)) + 1)
    step = (upper - lower) / n
    return [lower + n*step for n in range(n+1)], om


def show_folium(mymap,
                geodf,
                polygon_maker,
                description,
                value_unit,
                val_col,
                fill_color,
                threshold=0,
                bins_no=5,
                crs='EPSG:4326'):

    polygons = polygon_maker(geodf)
    mygdf = gpd.GeoDataFrame({
        'geometry': polygons, 
        val_col: geodf[val_col],
        'latitude': geodf.latitude,
        'longitude': geodf.longitude
    }, geometry='geometry', crs=crs)
    bins, om = make_bins(mygdf[val_col], bins_no)

    mygdf['evalue'] = mygdf[val_col] / math.pow(10, om)

    choropleth = folium.Choropleth(
        geo_data=mygdf[mygdf[val_col]>threshold],
        data=mygdf.evalue,
        key_on="feature.id",  # index of the geodataframe is transformed into `id` field
        fill_color=fill_color,
        fill_opacity=0.6,
        bins=bins,
        line_weight=0,
        line_color="black",
        nan_fill_opacity=0.0,
        legend_name=f"{val_col} [1.0e{om:+03d} {value_unit}]",  # title under the color scale
        name=f"{description}",  # name of thew layer, e.g. in the layer control
    ).add_to(mymap)

    cols_xy = [c for c in ["latitude", "longitude"] if c in mygdf]
    labels_xy = ["latitude:", "longitude:"] if cols_xy else []
    # add tooltip to appear, when pointing at a hectar
    tooltip = folium.GeoJsonTooltip(
        # column names with values to be displayed
        fields=["evalue"] + cols_xy,
        # text to be shown explaining each value
        aliases=[f"{val_col} [1.0e{om:+03d}]"] + labels_xy,
        localize=True,
        max_width=800,
    )
    tooltip.add_to(choropleth.geojson)

    return mymap


def map_maker(geodf, zoom_start=None):
    (sy,sx) = center(geodf)
    zoom_start = zoom_start or estimate_zoom_level(geodf)
    mymap = folium.Map(location=(sy,sx), zoom_start=zoom_start)
    return mymap


def exposures_folium(
        exposures,
        polygon_maker,
        add_to_map=None,
        val_col='value',
        threshold=0,
        bins_no=5,
        fill_color="PuBu",
        zoom_start=None):
    """Creates a folium (leaflet.js) map from the exposures object

    Parameters
    ----------
    exposures : climada.entity.Exposures
    polygon_maker : function(geopandas.GeoDataFrame) -> geopandas.GeoDataFrame
        applied to exposures.gdf
        supposed to make a copy and add a geometry column with Polygons around lat/lon
    add_to_map : folium.Map, optional
        the map where this layer is being integrated, 
        if None a new map is created based on the exposures' coordinates
    val_col : str, optional
        the column name of the exposures value of interest, default 'value'
    threshold : int, optional
        below this level values are not showed in the map, default 0
    bins_no : int, optional
        number of categories/bins to be shown, default 5
    fill_color : str, optional
        must be a ColorBrewer's name (see https://colorbrewer2.org)
        default: 'PuBu', Purple -> Blue
    zoom_start : int, optional
        starting zoom level, if None, the zoom level is guessed from the latitudinal range

    Returns
    -------
    folium.folium.Map
        A folium Map with the exposures values shading a polygon around their lat/lon coordinates
    """
    return show_folium(mymap=add_to_map or map_maker(exposures.gdf, zoom_start=zoom_start),
                geodf=exposures.gdf,
                polygon_maker=polygon_maker,
                description=exposures.description or val_col,
                value_unit=exposures.value_unit,
                val_col=val_col,
                fill_color=fill_color,
                threshold=threshold,
                bins_no=bins_no)


def impact_folium(
        impact,
        polygon_maker,
        add_to_map=None,
        label=None,
        threshold=0,
        bins_no=5,
        fill_color="YlOrBr",
        zoom_start=None):
    """Creates a folium (leaflet.js) map from the exposures object

    Parameters
    ----------
    impact : climada.engine.Impact
    polygon_maker : function(geopandas.GeoDataFrame) -> geopandas.GeoDataFrame
        applied to exposures.gdf
        supposed to make a copy and add a geometry column with Polygons around lat/lon
    add_to_map : folium.Map, optional
        the map where this layer is being integrated, 
        if None a new map is created based on the exposures' coordinates
    label : str, optional
        impact label, will be shown in the LayerControl of the folium.Map
        default is 'Expected annual impact for hazard type [HAZ_TYPE]'
    threshold : int, optional
        below this level values are not showed in the map, default 0
    bins_no : int, optional
        number of categories/bins to be shown, default 5
    fill_color : str, optional
        must be a ColorBrewer's name (see https://colorbrewer2.org)
        default 'YlOrBr', Yellow -> Orange -> Brown
    zoom_start : int, optional
        starting zoom level, if None, the zoom level is guessed from the latitudinal range

    Returns
    -------
    folium.folium.Map
        A folium Map with the exposures values shading a polygon around their lat/lon coordinates
    """
    imp_df = pd.DataFrame({
        'longitude': impact.coord_exp[:,1],
        'latitude': impact.coord_exp[:,0],
        'eai_exp': impact.eai_exp})
    label = label or f"Expected annual impact for hazard type {impact.haz_type}"
    return show_folium(mymap=add_to_map or map_maker(imp_df, zoom_start=zoom_start),
                geodf=imp_df,
                polygon_maker=polygon_maker,
                description=label,
                value_unit=impact.unit,
                val_col='eai_exp',
                fill_color=fill_color,
                threshold=threshold,
                bins_no=bins_no)


def add_layer_control(mymap):
    # https://github.com/python-visualization/folium/issues/816
    # "Both folium and Leaflet require the LayerControl to be added last: otherwise it doesn't know what layers there are."
    folium.LayerControl().add_to(mymap)
    return mymap
