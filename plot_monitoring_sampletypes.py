#!/usr/bin/env python

"""
Author: Lori Garzio on 7/15/2020
Last modified: 7/22/2020
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

    # add MARCO cruise data
    # marco = 'AcidificationMonitoringMidA_Ver201902.csv'
    # marco_df = pd.read_csv(os.path.join(file_dir, marco))
    #
    # df_ct1 = marco_df[(marco_df['Collection_Type'] == 'Cruise') |
    #                   (marco_df['Collection_Type'] == 'Continuous Monitoring') |
    #                   (marco_df['Collection_Type'] == 'Ongoing Fixed Station')]
    #
    # # remove Maryland data points
    # df_ct = df_ct1[(df_ct1['X'] > datacoords['minlon']) & (df_ct1['X'] < datacoords['maxlon']) &
    #                (df_ct1['Y'] > datacoords['minlat']) & (df_ct1['Y'] < datacoords['maxlat']) &
    #                (df_ct1['Entity'] != 'MDDNR - ConMon') & (df_ct1['Entity'] != 'MDDNR - DFlow') &
    #                (df_ct1['Entity'] != 'MDDNR - MDcoastalBays')]
    #
    # add_col = []
    # oa_samples = ['pH', 'DIC', 'TA', 'pCO2']
    # for index, row in df_ct.iterrows():
    #     n = 0
    #     for oas in oa_samples:
    #         if row[oas].lower() == 'yes':
    #             n += 1
    #
    #     if n > 1:
    #         add_col.append('2+ OA samples')
    #         # df_ct.loc[index, 'oa_sampling'] = '2+ OA samples'
    #     else:
    #         for oas in oa_samples:
    #             if row[oas].lower() == 'yes':
    #                 add_col.append('{} only'.format(oas))
    #                 # df_ct.loc[index, 'oa_sampling'] = '{} only'.format(oas)
    #
    # df_ct['oa_sampling'] = add_col
    # df_ct.to_csv('/Users/lgarzio/Documents/rucool/Saba/NJDEP_OA2020/df_ct.csv')
    df_ct = pd.read_csv('/Users/lgarzio/Documents/rucool/Saba/NJDEP_OA2020/df_ct.csv')
    sampling = np.unique(df_ct['oa_sampling'])
    sampling = ['pH only', 'pCO2 only', '2+ OA samples']
    colors = ['black', 'cyan', 'gold']
    #rc('text', usetex=True)  # activate latex text rendering

    for i, s in enumerate(sampling):
        idf = df_ct[df_ct['oa_sampling'] == s]
        if s == '2+ OA samples':
            lab = '2+ parameters'
        elif s == 'pCO2 only':
            lab = ''.join(('$\it{p}$', r'$\rmCO_2$'))
        else:
            lab = s
        ax.scatter(idf['X'], idf['Y'], color=colors[i], s=16, edgecolor='black', linewidths=.2,
                   label=lab, zorder=10)

    # add EcoMon data
    emf = 'EcoMon_OA_stations_mod.csv'
    df = pd.read_csv(os.path.join(file_dir, emf))

    emlat = np.array(df['lat_deg']) + np.array(df['lat_min']) / 60
    emlon = -(np.array(df['lon_deg']) + np.array(df['lon_min']) / 60)

    # only plot data points within data coordinates
    lon_ind = np.logical_and(emlon > datacoords['minlon'], emlon < datacoords['maxlon'])
    lat_ind = np.logical_and(emlat > datacoords['minlat'], emlat < datacoords['maxlat'])

    coordind = np.where(np.logical_and(lon_ind, lat_ind))

    ax.scatter(emlon[coordind], emlat[coordind], color='gold', s=16, edgecolor='black', linewidths=.2, zorder=10)

    # add glider tracks
    gl = 'gliders.csv'
    df = pd.read_csv(os.path.join(file_dir, gl))
    gliders = ['a', 'b']

    for i, g in enumerate(gliders):
        dfg = df[df['Name'] == g]
        gllat = np.array(dfg['lat'])
        gllon = np.array(dfg['lon'])
        ax.plot(gllon, gllat, c='black', lw=2, zorder=10)

    plt.legend(fontsize=10, loc='upper left', facecolor='white', framealpha=1)

    plt.savefig(os.path.join(file_dir, 'OA_monitoring_map_details.png'), dpi=200)
    plt.close()


if __name__ == '__main__':
    fpath = '/Users/lgarzio/Documents/rucool/Saba/NJDEP_OA2020'
    main(fpath)
