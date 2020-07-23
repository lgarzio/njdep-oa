#!/usr/bin/env python

"""
Author: Lori Garzio on 7/22/2020
Last modified: 7/23/2020
Buoy locations from https://www.ndbc.noaa.gov/
Shellfish hatchery list from https://ecsga.org/growers-resources/
Glider data from http://slocum-data.marine.rutgers.edu/erddap/search/advanced.html?page=1&itemsPerPage=1000&searchFor=ru28-20190717T1522
"""

import numpy as np
import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console
plt.rcParams.update({'font.size': 13})


def add_map_features(ax, axes_limits, land_color=None):
    """
    Adds latitude and longitude gridlines and labels, coastlines, and statelines to a cartopy map object
    :param ax: plotting axis object
    :param axes_limits: list of axis limits [min lon, max lon, min lat, max lat]
    :param land_color: optional color for land
    """
    #gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='dotted', x_inline=False)
    gl = ax.gridlines(draw_labels=True, linewidth=0)
    gl.top_labels = False
    gl.right_labels = False
    #gl.rotate_labels = False
    gl.xpadding = 12
    gl.ypadding = 12
    ax.set_extent(axes_limits)

    if land_color:
        lc = land_color
    else:
        lc = cfeature.COLORS['land_alt1']

    land = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='black', facecolor=lc)
    ax.add_feature(land, zorder=1)

    state_lines = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(state_lines, edgecolor='black')

    feature = cfeature.NaturalEarthFeature(
        name='coastline', category='physical',
        scale='10m',
        edgecolor='black', facecolor='none')
    ax.add_feature(feature, zorder=8)

    rivers = cfeature.NaturalEarthFeature(
        category='physical', name='rivers_lake_centerlines',
        scale='10m', facecolor='none', edgecolor='lightblue')
    ax.add_feature(rivers, linewidth=1)


def set_map(data):
    """
    Set up the map and projection
    :param data: data from the netcdf file to be plotted, including latitude and longitude coordinates
    :returns fig, ax objects
    :returns dlat: latitude data
    returns dlon: longitude data
    """
    #lccproj = ccrs.LambertConformal(central_longitude=-74.5, central_latitude=38.8)
    #fig, ax = plt.subplots(figsize=(8, 9), subplot_kw=dict(projection=lccproj))
    fig, ax = plt.subplots(figsize=(8, 9), subplot_kw=dict(projection=ccrs.PlateCarree()))

    dlat = data['XLAT'].values
    dlon = data['XLONG'].values

    return fig, ax, dlat, dlon


def main(file_dir):
    pltcoords = dict(minlon=-76, maxlon=-71.4, minlat=38, maxlat=41.8)
    datacoords = dict(minlon=-75.5, maxlon=-71.6, minlat=38.1, maxlat=41.8)

    # set up figure
    fig, ax = plt.subplots(figsize=(8, 9), subplot_kw=dict(projection=ccrs.PlateCarree()))
    axlims = [pltcoords['minlon'], pltcoords['maxlon'], pltcoords['minlat'], pltcoords['maxlat']]
    add_map_features(ax, axlims, 'darkgray')
    #add_map_features(ax, axlims)

    # add bathymetry
    gf = os.path.join(file_dir, 'GMRTv3_7_20200716topo-mask.grd')
    grid_file = xr.open_dataset(gf)
    gf_lon = grid_file['lon']
    gf_lat = grid_file['lat']
    gf_lon_ind = np.logical_and(gf_lon > datacoords['minlon']-1, gf_lon < datacoords['maxlon']+1)
    gf_lat_ind = np.logical_and(gf_lat > datacoords['minlat']-1, gf_lat < datacoords['maxlat']+1)
    bathy = grid_file['altitude'][gf_lat_ind, gf_lon_ind].values

    land_mask = np.logical_and(bathy >= -1, bathy >= -1)
    bathy[land_mask] = np.nan  # turn values over land to nans
    depth = abs(bathy)
    ln_depth = np.log(depth)
    levs = [1, 25, 50, 100, 500, 1000, 2000, 3000, 4000]
    ln_levs = np.log(levs)
    #ax.contour(gf_lon[gf_lon_ind], gf_lat[gf_lat_ind], depth, levs, colors='black', linestyles='solid', linewidths=.75)
    cs = ax.contourf(gf_lon[gf_lon_ind], gf_lat[gf_lat_ind], ln_depth, ln_levs, cmap='Blues', alpha=.75)

    # add NJ shellfish hatcheries
    hatch_names = ['4C', 'Bayfarm', 'BillAvery', 'ClamDaddy', 'IntracoastalAquaculture', 'NauticalNuggets',
                   'NJAquacultureInnovation', 'RutgersCapeShore']
    hatch_lons = [-74.659, -74.330, -74.454, -74.365, -74.340, -74.486, -74.942, -75.031]
    hatch_lats = [39.195, 39.544, 39.469, 39.410, 39.603, 39.521, 38.969, 39.234]

    ax.scatter(hatch_lons, hatch_lats, color='cyan', marker='s', s=50, edgecolor='black', label='Hatcheries',
               zorder=10)

    # add NOAA buoy locations
    noaa_names = ['44065', '44025', '44066', '44091', '44009']
    noaa_lons = [-73.703, -73.164, -72.644, -73.769, -74.702]
    noaa_lats = [40.369, 40.251, 39.618, 39.778, 38.457]

    ax.scatter(noaa_lons, noaa_lats, color='yellow', marker='D', s=50, edgecolor='black', label='NOAA buoys', zorder=10)

    # for i, n in enumerate(noaa_names):
    #     ax.text(noaa_lons[i] + .1, noaa_lats[i] + .02, n)

    # add NERRS buoy locations
    nerrs = ['JCTN4', 'JCQN4']
    nerrs_lons = [-74.338, -74.461]
    nerrs_lats = [39.508, 39.548]

    ax.scatter(nerrs_lons, nerrs_lats, color='orange', marker='.', s=175, edgecolor='black', label='JCNERR buoys', zorder=10)

    # for i, n in enumerate(nerrs):
    #     if n == 'JCTN4':
    #         ax.text(nerrs_lons[i] + .15, nerrs_lats[i] - .02, n)
    #     elif n == 'JCQN4':
    #         ax.text(nerrs_lons[i] - .2, nerrs_lats[i] + .05, n)

    # add glider
    gf = '/Users/lgarzio/Documents/rucool/Saba/NJDEP_OA2020/ru28-20190717T1522-profile-sci-rt_cbd2_aed9_b2b9.nc'
    ds = xr.open_dataset(gf, mask_and_scale=False)
    gllat = ds['latitude'].values
    gllon = ds['longitude'].values

    ax.plot(gllon, gllat, color='magenta', label='Glider')

    plt.legend(fontsize=12, loc='upper left', facecolor='white', framealpha=1)

    plt.savefig(os.path.join(file_dir, 'OA_proposed_sites_draft.png'), dpi=200)
    plt.close()


if __name__ == '__main__':
    fpath = '/Users/lgarzio/Documents/rucool/Saba/NJDEP_OA2020'
    main(fpath)
