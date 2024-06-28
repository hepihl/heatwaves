'''
A set of functions that perform routine calculations/operations for basic atmospheric/ocean data 
visualization.

6/10/2024
Author: Haakon Pihlaja
'''

from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy
import xarray as xr
import scipy
import geopy.distance

def anomalies(data, gb_time="time.month", rs_time="1MS", monthly=False):
    """
    Calculates anomalies in given variable using 1991-2020 standard.

    :data (DataArray): unmanipulated environmental data
    :time (str): string specifying the window used to resample data
    :monthly (bool): boolean value specifying whether input data is already a monthly mean
    :return (DataArray): anomalies
    """
    clim = data.sel(time=slice("1991-01-01", "2020-12-31"))
    
    if monthly==False:
        anom = data.resample(time=rs_time).mean(dim="time").groupby(gb_time) - clim.groupby(gb_time).mean(dim="time")
    else:
        anom = data.groupby("time.month") - clim.groupby("time.month").mean(dim="time")
    
    return anom

def plt_temp(data, ax, *args, **kwargs):
    """
    Creates plot of sea surface temperatures in the North Atlantic (0 to 60 N, 0 to 80 W). 

    :data (DataArray): SST data
    :ax (matplotlib axes): axis onto which the SSTs will be plotted
    """
    if ax.has_data():
        temp = data.plot(ax=ax, *args, **kwargs)
        return temp, ax
    else:
        central_lon = -40
        central_lat = 30
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAKES, edgecolor='black', facecolor="none")
        ax.set_extent([-80, 0, 0, 60])
        gl = ax.gridlines(alpha = 0.5, draw_labels=True, linestyle = "--", x_inline=False, y_inline=False)
        gl.right_labels = False
        gl.top_labels = False
        temp = data.plot(ax=ax, *args, **kwargs)
        return temp, ax
    
def plt_wind(data, ax, *args, **kwargs):
    """
    Creates quiver plot of 10 m wind field over the North Atlantic (0 to 60 N, 0 to 80 W). 

    :data (DataArray): wind data (you'll want to downsample this using xarray's .coarsen() method)
    :ax (matplotlib axes): axis onto which the winds will be plotted
    """
    if ax.has_data():
        data.plot.quiver(ax=ax, *args, **kwargs)
        return ax
    else:
        central_lon = -40
        central_lat = 30
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAKES, edgecolor='black', facecolor="none")
        ax.set_extent([-80, 0, 0, 60])
        gl = ax.gridlines(alpha = 0.5, draw_labels=True, linestyle = "--", x_inline=False, y_inline=False)
        gl.right_labels = False
        gl.top_labels = False
        data.plot.quiver(ax=ax, *args, **kwargs)
        return ax

def plt_pressure(data, ax, *args, **kwargs):
    """
    Creates contour plot of sea level pressure over the North Atlantic (0 to 60 N, 0 to 80 W). 

    :data (DataArray): pressure data
    :ax (matplotlib axes): axis onto which the presure contours will be plotted
    """
    if ax.has_data():
        pressure = data.plot.contourf(ax=ax, *args, **kwargs)
        return pressure, ax
    else:
        central_lon = -40
        central_lat = 30
        ax.coastlines()
        ax.add_feature(cartopy.feature.LAKES, edgecolor='black', facecolor="none")
        ax.set_extent([-80, 0, 0, 60])
        gl = ax.gridlines(alpha = 0.5, draw_labels=True, linestyle = "--", x_inline=False, y_inline=False)
        gl.right_labels = False
        gl.top_labels = False
        pressure = data.plot.contourf(ax=ax, *args, **kwargs)
        return pressure, ax

def monthly_plt(data, axs=[], quiver=False, contour=False):
    """
    Creates a 12-month plot of given data. 
    6/10/2024 - I wouldn't recommend using this... it was a first-pass attempt at streamlining 
                my code, but it't not very robust. See atm_sst.ipynb for a better setup.

    :data (DataArray): environmental data, generally anomalies
    :return (plt.figure, plt.axs): figure and axes
    """

    data = data.groupby("time.month").mean(dim="time")
    central_lon = -40
    central_lat = 30
    coloring = plt.get_cmap("seismic")
    months = {1:"January", 2:"February", 3:"March", 4:"April", 5:"May", 6:"June", 7:"July", 8:"August", 9:"September", 10:"October", 11:"November", 12:"December"}
    
    if len(axs)==0:
        fig, axs = plt.subplots(3, 4, figsize=(7.2, 5), layout = 'compressed', subplot_kw={"projection":ccrs.LambertConformal(central_lon, central_lat)})
        axs = axs.flatten()
    
        for i,month in enumerate(range(1,13)):
            anom_month = data.sel(month=month)
            anom_month["month"] = months[i+1]
            axs[i].coastlines()
            axs[i].add_feature(cartopy.feature.LAKES, edgecolor='black', facecolor="none")
            axs[i].set_extent([-80, 0, 0, 60])
            gl = axs[i].gridlines(alpha = 0.5, draw_labels=True, linestyle = "--", x_inline=False, y_inline=False)
            gl.right_labels = False
            gl.top_labels = False
            gl.xlabel_style = {"size":6}
            gl.ylabel_style = {"size":6}
            anom_plt = anom_month.plot(ax=axs[i], cmap=coloring, transform=ccrs.PlateCarree(), vmin=-5, vmax=5, add_colorbar=False)
            axs[i].set_title(anom_month["month"].values, fontsize=10)
    
        cbar = plt.colorbar(anom_plt, ax=axs, orientation="vertical", extend="both", aspect=30, shrink=0.6)
        cbar.ax.set_ylabel('SST Anomaly ($^\circ$C)', fontsize=10)
        return fig, axs
        
        
    elif quiver == True:
        for i,month in enumerate(range(1,13)):
            anom_month = data.sel(month=month)
            anom_plot = anom_month.plot.quiver(ax=axs[i], x="lon", y="lat", u="u10", v="v10", transform=ccrs.PlateCarree(), add_guide=False)
            axs[i].set_title(anom_month["month"].values, fontsize=10)

        return anom_plot
            
    elif contour == True:
        for i,month in enumerate(range(1,13)):
            anom_month = data.sel(month=month)
            anom_plt = anom_month.plot(ax=axs[i], cmap=coloring, transform=ccrs.PlateCarree(), vmin=-5, vmax=5, add_colorbar=False)
            axs[i].set_title(anom_month["month"].values, fontsize=10)

        return anom_plot

def ohc_depth(temp, depth, cw=4180, p0= 1000):
    """
    Calculates depth-integrated ocean heat content 

    :temp (DataArray): potential temperature, must contain dimensions "time", "depth", "lat", and "lon"
    :depth (str or int): depth over which to integrate, use "full" to integrate over full depth
    :cw (int): specific heat capacity of seawater (J kg^-1 K^-1), default value taken from Marshall and Plumb
    :p0 (int): reference density (kg m^-3)
    :return (DataArrays): depth-intergrated ocean heat content (J/m^2)
    """
    # Need to fill missing values, e.g. land or sea floor, with 0
    # Otherwise NaNs will propagate through integration
    temp = temp.fillna(0) 
    
    if depth == "full":           
        integral = xr.apply_ufunc(scipy.integrate.trapezoid, temp, kwargs={"x":temp["depth"]},input_core_dims = [["depth"]])
        ohc = cw * p0 * integral
    else:
        temp = temp.sel(depth=slice(0, depth))
        integral = xr.apply_ufunc(scipy.integrate.trapezoid, temp, kwargs={"x":temp["depth"]},input_core_dims = [["depth"]])
        ohc = cw * p0 * integral

    return ohc

def ohc_horiz(data, to_int):
    """
    Calculates horizontally-integrated ocean heat content. Assumes data has already been depth-integrated.

    :data (DataArray): depth-integrated ocean heat content, must contain dimensions "time", "lat", and "lon"
    :to_int (list): list of dimensions over which to integrate, should be one of ['lat'], ['lon'], ['lat','lon']
    :return (DataArrays): horizontally-integrated ocean heat content (J)
    """
    if "lat" not in to_int:
        q_bar_x = xr.DataArray()
        xs = []

        # Calculate distance between lines of longitude at each latitude there is data for
        for lat in data["lat"].values:
            pt1 = [lat, data["lon"].values[0]]
            pt2 = [lat, data["lon"].values[1]]
            dist = geopy.distance.distance(pt1, pt2)
            xs.append(dist.m) 
        
        # Integrate along lines of latitude, taking dx to be the distances calculated in the above loop
        for i in range(len(data['lat'].values)):
            zonal_int = xr.apply_ufunc(scipy.integrate.trapezoid, data.isel(lat=i), kwargs={"dx":xs[i]},input_core_dims = [["lon"]])
            if i == 0 :
                q_bar_x = zonal_int
            else:
                q_bar_x = xr.concat([q_bar_x, zonal_int], dim="lat")

        return q_bar_x

    elif "lon" not in to_int:
        # Integrate along lines of longitude, taking dx to be a constant 111 km
        q_bar_y = xr.apply_ufunc(scipy.integrate.trapezoid, data, kwargs={"dx":111000},input_core_dims = [["lat"]])
        return q_bar_y

    else:
        q_bar_x = xr.DataArray()
        xs = []

        # Calculate distance between lines of longitude at each latitude there is data for
        for lat in data["lat"].values:
            pt1 = [lat, data["lon"].values[0]]
            pt2 = [lat, data["lon"].values[1]]
            dist = geopy.distance.distance(pt1, pt2)
            xs.append(dist.m) 
        
        # Integrate along lines of longitude, taking dx to be the distances calculated in the above loop
        for i in range(len(data['lat'].values)):
            zonal_int = xr.apply_ufunc(scipy.integrate.trapezoid, data.isel(lat=i), kwargs={"dx":xs[i]},input_core_dims = [["lon"]])
            if i == 0 :
                q_bar_x = zonal_int
            else:
                q_bar_x = xr.concat([q_bar_x, zonal_int], dim="lat")

        # Integrate along lines of latitude
        q_bar = xr.apply_ufunc(scipy.integrate.trapezoid, q_bar_x, kwargs={"dx":111000},input_core_dims = [["lat"]])
        return q_bar