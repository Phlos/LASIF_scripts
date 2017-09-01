#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plotting a map that gives the misfit per station

:copyright:
  Nienke Blom (n.a.blom@uu.nl), 2016, 2017
  
"""
import os
from lasif.components.project import Project
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np
from matplotlib import cm

def set_axis_in_pi_units(ax, axtype):

  # determine current maximum of axis
  if axtype == 'cbar_x' or axtype == 'cbar_y':
    lims = ax.get_clim()
  elif axtype == 'x':
    lims = ax.get_xlim()
  elif axtype == 'y':
    lims = ax.get_ylim()
  else:
    print 'unrecognised axis type: {}'.format(axtype)
    return

  axlength = lims[1] - lims[0]
  # determine how many pi (sub)units fit in there
  tickbase = np.pi; divisor=2; prev_divisor=1
  while axlength / tickbase < 3:
    tickbase /= divisor
    tickbase *= prev_divisor
#    print 'new tickbase for {:.2f}:   {:.3f}*pi, \
#           how many of them? {:.1f}. \
#           Divisor is {}, before was {}'.format(axlength, 
#                                                tickbase/np.pi, 
#                                                axlength/tickbase, 
#                                                divisor, prev_divisor)
    prev_divisor = divisor
    if divisor == 2:
      divisor = 4
    elif divisor == 4:
      divisor = 5
    elif divisor == 5:
      divisor = 10
    elif divisor == 10:
      divisor = 2; prev_divisor = 1

  # set the ticks such 
  def piify_labels(a):
    for _ii, x in enumerate(a):
      x =  float(x) / np.pi
      x =  '{:g}'.format(x)+u" \u03c0"
      a[_ii] = x
    return a
    
  def piify_labels_cbar(a):
    for x in a:
      tx = str(x.get_text().replace(u"\u2212", '-'))
      tx = float(tx) / np.pi
      tx = '{:.0e}'.format(tx)
      tx = float(tx)
      tx = '{:g}'.format(tx)+u" \u03c0"
      x.set_text(tx)
    return a
  
  if axtype == 'x':
    #ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    ax.xaxis.set_major_locator(tck.MultipleLocator(base=tickbase))
    a = piify_labels(ax.get_xticks().tolist())
    ax.set_xticklabels(a)
  elif axtype == 'y':
    #ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    ax.yaxis.set_major_locator(tck.MultipleLocator(base=tickbase))
    a = piify_labels(ax.get_yticks().tolist())
    ax.set_yticklabels(a)
  elif axtype == 'cbar_x':
    muloc = tck.MultipleLocator(base=tickbase)
    ax.set_ticks(muloc.tick_values(lims[0], lims[1]))
    a = piify_labels_cbar(ax.ax.get_xticklabels())
    ax.ax.set_xticklabels(a)
  elif axtype == 'cbar_y':
    muloc = tck.MultipleLocator(base=tickbase)
    ax.set_ticks(muloc.tick_values(lims[0], lims[1]))
    a = piify_labels_cbar(ax.ax.get_yticklabels())
    ax.ax.set_yticklabels(a)
    

def plot_title_legend(event_info, nwinsta, nsta, itername):
  event_name = event_info['event_name']
  where = event_info['region']
  when = str(event_info['origin_time'].strftime('%Y-%m-%d %H:%M'))
  mag   = str(event_info['magnitude'])
  magtype = event_info['magnitude_type']
  titel='Event '+event_name+' ('+when+') - '+ \
	mag+' '+magtype+' - '+ \
	str(nwinsta)+' stations with windows ('+ \
	str(round(100*nwinsta/nsta))+'% of '+str(nsta)+') - '+itername
  plt.title(titel)
  plt.legend()

def plot_coloured_hist(ax, plotted_shite, colrange, colmap):
  
  '''
  plot histograms that are coloured by colmap in the range specified by colrange
  '''
  
  minn = colrange[0]; maks = colrange[1]
  n, bins, patches = ax.hist(np.clip(plotted_shite, minn, maks), 50, range=colrange, facecolor='green', alpha=0.75)
  ax.set_xlim(minn, maks)
  
  # colour-code histogram
  bin_centers = 0.5 * (bins[:-1] + bins[1:])
  # scale values to interval [0,1]
  col = bin_centers - min(bin_centers)
  col /= max(col)
  for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', colmap(c))
    
def plot_misfits_map(comm, iteration_name, event_name, plottype, 
                     lngs, lats, lngs_nan, lats_nan, 
                     plotted_shite, colmap, norm, colrange, collabel,
                     plot_windowless, plot_beachball, plot_histogram,
                     verbose=False):
  
  '''
  Plots misfits (unspecified, but gathered in plotted_shite)
  on a map, together with a histogram (or not)
  '''
  
  event_info = comm.events.get(event_name)
  
  # prepare figure: one or two subplots
  if plot_histogram:
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(12.25,12.75),
                                   gridspec_kw={'height_ratios':[5,1]}#,
                                                #'width_ratios':[1,1]}
                                  )
    #fig.set_size_inches(10., 18., forward=True)
  else:
    fig, ax1 = plt.subplots(1, figsize=(18.,13.5))
    #fig.set_size_inches(18.0, 13.5, forward=True)




  #----------------------
  #--- the map
  #----------------------
  
  # make basic LASIF style map (w/ phyical domain plotted)
  
  if verbose:
    print 'Making map object'
  map_object = comm.project.domain.plot(ax=ax1)

  # convert lat/lon to map-projected values
  if verbose:
    'Converting lat/lon'
  x, y = map_object(lngs, lats)
  xnan,ynan = map_object(lngs_nan, lats_nan)

  #- actually plot the misfits

  if verbose:
    print 'Plotting...'
  if plot_windowless:
    windowless = map_object.scatter(xnan, ynan, color='black', 
                                    s=50, marker="v", alpha=0.2,
                                    label="stations without windows")

  # plot the actual misfits.
  misfitplot = map_object.scatter(x, y, c=plotted_shite, 
                                  norm=norm, cmap=colmap, 
                                  s=50, marker="o",
                                  alpha=0.7, zorder=5,
                                  label=plottype)

  # plot the event.
  if plot_beachball == True:
    from lasif import visualization
    visualization.plot_events(events=[event_info], map_object=map_object)

  # add colour bar.
  cbar = map_object.colorbar(misfitplot,location='bottom',pad="5%")
  if plottype is not 'misfits':
    set_axis_in_pi_units(cbar, 'cbar_x')
  cbar.set_label(collabel)

  # add title
  nsta = len(lngs) + len(lngs_nan)
  nwinsta = len(lngs)
  plot_title_legend(event_info, nwinsta, nsta, 'iteration '+iteration_name)
  
  
  
  #----------------------
  #--- the histogram
  #----------------------

  if plot_histogram:
    plot_coloured_hist(ax=ax2, plotted_shite=plotted_shite, colrange=colrange, colmap=colmap)
    if plottype == 'misfits':
      colr = '[{:.2g}, {:.2g}]'.format(colrange[0], colrange[1])
    else:
      colr = u"[{:.2g} \u03c0, {:.2g} \u03c0]".format(colrange[0]/np.pi, colrange[1]/np.pi)
    ax2.set_xlabel('station averaged '+plottype.replace('_', ' ')+' (values beyond '+colr+' binned in last bins!)')
    if plottype is not 'misfits':
      set_axis_in_pi_units(ax=ax2, axtype='x')

  return fig




def plot_misfit_map_avg(Lasif_path='./EMed_full', 
                        iteration_name='10_new', 
                        event_name='event_01',
                        plottype = "phaseshifts",   # or 'misfits'
                        plot_beachball=True, 
                        plot_windowless=True,
                        plot_histogram=True,
                        colour_max=np.nan,
                        showdontsave=True,
                        savepath='./plots/',
		                    output=True, verbose=True, veryverbose=False):
  '''
  
  Plots total misfits (or some other metric) per station of 
  a specific event at a specific iteration on a map, where 
  the misfit value is shown as a colour scale. The stations 
  without windows can be plotted in transparent triangles, 
  and additionally a histogram can be included that shows 
  the distribution of misfits / phase shifts. 
  
  plottype = 'misfits':      the TOTAL misfit per station is plotted
  plottype = 'phasediff_integral': the AVERAGE phaseshift integral is plotted
  plottype = 'mean_delay':   the mean delay for the whole window, averaged over all windows
  plottype = 'max_delay':    the max delay (+ or -) for the whole window, averaged over all windows
  plottype = 'median_delay': the median (middle value) delay, averaged over all windows
  
  INPUT:
  Lasif_path='./EMed_full'  # path to lasif project
  iteration_name = "10_new" # Lasif iteration name for which misfits should be plotted
  event_name = "event_01"   # Lasif event name of event to be plotted
  plot_beachball = True     # plot the beachball of the earthquake
  plot_windowless = True    # plot stations without windows or not
  plottype = "phaseshifts"  # can also be "misfits"
  colour_max=[any number]   # determines the colour scale of the misfits
  showdontsave=True         # will show the picture instead of saving it
  savepath='./plots/'       # location where the plots are saved
  '''
  
  valid_plottypes=['phasediff_integral', 'misfits', 'mean_delay', 'median_delay', 'max_delay']
  
  if plottype not in valid_plottypes:
    print 'Warning: plottype must be one of '+', '.join(valid_plottypes)
    return

  # import the project communicator
  comm = Project(Lasif_path).comm

  # obtain information about station locations
  stations = comm.query.get_all_stations_for_event(event_name)
  nsta = len(stations)

  # obtain information about the windows
  wm = comm.windows.get(event_name, iteration_name)
  # station list for windows
  window_channels = wm.list()

  if len(window_channels) == 0:
    print "Event "+event_name+" contained no windows in iteration "+iteration+" in project "+Lasif_path
    
  #- obtain station info and misfit info
  #  (I copied/adapted this from plot_stations_for_event)
  lngs = []
  lats = []
  lngs_nan = []
  lats_nan = []
  station_ids = []
  misfits = []
#  phaseshifts  = [] #  can I immediately make ndarrays out of these? 
                    #+ Saves some converting later on when the faffing 
                    #+ with leaving list elements out is relevant (the isnan stuff)).
  visualisation_metric = []
  
  # initialise the saving of extremal values
  extreme_pos = []
  extreme_neg = []

  if veryverbose:
    print 'Calculating misfits etc'
  #- for every station for that particular event, obtain misfits
  for stationID, value in stations.iteritems():  # Loop as dicts are unordered.
      
    # get windows for station
    wins = wm.get_windows_for_station(stationID)
    
    # get window info
    # (first determine whether there are windows at all (also fuck ugly))
    sz = np.size(wins) # sz will always be max 3: for each channel Z/N/E one
    if sz > 0:

      # get the misfits of all windows
      misfitsum = 0; phshiftsum = 0; viz_metric_sum = 0
      nwins = 0
      for ii in np.arange(sz):
        channelID = wins[ii].channel_id
        
        for jj in np.arange(len(wins[ii].windows)):
          nwins += 1
#          print 'channel {}: {}, window {}'.format(ii,channelID, jj)
          
          # if plotting misfit...
          if plottype == 'misfits':
            misfit_val = wins[ii].windows[jj].misfit_value # non-normalised
            misfit_wt  = wins[ii].windows[jj].weight # window weight
            misfit_weighted = misfit_val * misfit_wt 
            misfitsum += misfit_weighted
          
          # if plotting something else (phaseshift_integral, mean_delay, median_delay)
          else:
            if plottype in wins[ii].windows[jj].misfit_details:
              viz_metric = float(wins[ii].windows[jj].misfit_details[plottype])
            else:
              print 'Warning: no '+plottype+' saved for channel '+channelID+', window {} !'.format(jj)
              viz_metric = 0;
            viz_metric_sum += viz_metric

      # calculate _average_ of metric over all windows
      viz_metric_avg = viz_metric_sum/nwins

      # add misfit/phaseshift to array
      misfits.append(misfitsum)
      visualisation_metric.append(viz_metric_avg)
      
      # add lonlats to arrays
      lngs.append(value["longitude"])
      lats.append(value["latitude"])
      station_ids.append(stationID)

      if not np.isnan(colour_max):
        if plottype=='misfits':
          if misfitsum > colour_max:
            extreme_pos.append(stationID)
          elif misfitsum < 0:
            extreme_neg.append(stationID)
        else:
          if viz_metric_avg > colour_max:
            extreme_pos.append(stationID)
          elif viz_metric_avg < -colour_max:
            extreme_neg.append(stationID)
      
    else:
      lngs_nan.append(value["longitude"])
      lats_nan.append(value["latitude"])

  nwinsta = len(lngs)
  
  
  #- prepare the plotted shite:
  #   - put in a ndarray as otherwise indexing is bothersome
  #   - also make a color scale with 0 either as minimum or middle value
  import matplotlib.colors as colors
  from matplotlib import cm
  if plottype == "misfits":
    plotted_shite = np.asarray(misfits)
    if np.isnan(colour_max):
      maks = np.max(np.abs( plotted_shite[~np.isnan(plotted_shite)] ))
    else:
      maks = colour_max
    minn = 0
    norm = colors.Normalize(vmin=minn, vmax=maks)
    colmap = cm.Reds
    collabel = "total misfit per station (TF phase misfit, Fichtner et al. GJI 2008)"
    
  else:
    plotted_shite = np.asarray(visualisation_metric)
    if np.isnan(colour_max):
      maks = np.max(np.abs( plotted_shite[~np.isnan(plotted_shite)] ))
    else:
      maks = colour_max
    minn = -maks
    norm = colors.Normalize(vmin=minn, vmax=maks)
    colmap = cm.RdBu
    collabel = "station averaged {} (TF phase misfit, Fichtner et al. GJI 2008)".format(plottype.replace('_', ' '))
    
  print "colour maximum is "+str(maks)  
  colrange = [minn, maks]
  
  
  
  ###################################################
  ## PLOTTING
  ###################################################
  
  fig = plot_misfits_map(comm, iteration_name, event_name, plottype, 
                   lngs, lats, lngs_nan, lats_nan, 
                   plotted_shite, colmap, norm, colrange, collabel,
                   plot_windowless, plot_beachball, plot_histogram, verbose=veryverbose)

  # add total misfit value to histogram
  if plot_histogram and plottype=='misfits':
    #print 'blips'
    (ax1, ax2, ax3) = fig.get_axes()
    total_value = plotted_shite.sum()
    ax2.text(0.01, 0.85, 'Total misfit for event '+event_name+', iter '+iteration_name+': '+'{0:.3f}'.format(total_value), transform=ax2.transAxes)
  # save & show
  #plt.tight_layout()
  if showdontsave:
    plt.ion()
    plt.show()
    plt.ioff()
  else:
    savename = 'misfits.'+iteration_name+'.'+event_name+'.'+plottype+'.avg'
    if plot_windowless:
      savename += '.with-winless'
    if plot_histogram:
      savename += '.with-hist'
    savename += '.png'
    savepath = os.path.join(savepath, savename)
    #print 'saving to: '+savepath
    plt.savefig(savepath)
    print 'saved figure to '+savepath
    plt.close()
  
  extreme_neg = sorted(extreme_neg)
  extreme_pos = sorted(extreme_pos)
  extreme_val = sorted(extreme_neg + [s + '  +' for s in extreme_pos])
  
  if verbose:
    print 'stations with extremely strong '+plottype+' are:\n'+'\n'.join(extreme_val)
  
  if output:
    return extreme_neg, extreme_pos


def plot_misfit_change_map(Lasif_path='./EMed_full', 
                           iteration1='0_new', 
                           iteration2='9_4',
                           event_name='event_01',
                           plottype = "phaseshifts",   # or 'misfits'
                           plot_beachball=True, 
                           plot_windowless=True,
                           plot_histogram=False,
                           colour_max=np.nan,
                           showdontsave=True,
                           savepath='./plots/',
		                       output=False, verbose=True, veryverbose=False):
  '''
  
  Plots the change in misfit (or phaseshift) of a specific event 
  at a specific iteration on a map, where the misfit value is
  shown as a colour scale. The stations without windows can
  be plotted in transparent triangles, and additionally a 
  histogram can be included that shows the distribution of 
  misfits / phase shifts.
  
  plottype = 'misfits':      the TOTAL misfit per station is plotted
  plottype = 'phasediff_integral': the AVERAGE phaseshift integral is plotted
  plottype = 'mean_delay':   the mean delay for the whole window, averaged over all windows
  plottype = 'median_delay': the median (middle value) delay, averaged over all windows
  
  INPUT:
  Lasif_path='./EMed_full'  # path to lasif project
  iteration1 = "0_new"      # Lasif iteration name for first iteration to be compared
  iteration2 = "9_4"        # Lasif iteration name for second iteration to be compared
  event_name = "event_01"   # Lasif event name of event to be plotted
  plot_beachball = True     # plot the beachball of the earthquake
  plot_windowless = True    # plot stations without windows or not
  plottype = "phaseshifts"  # can also be "misfits"
  colour_max=[any number]   # determines the colour scale of the misfits
  showdontsave=True         # will show the picture instead of saving it
  savepath='./plots/'       # location where the plots are saved
  '''
  
  valid_plottypes=['phasediff_integral', 'misfits', 'mean_delay', 'median_delay', 'max_delay']

  if plottype not in valid_plottypes:
    print 'Warning: plottype must be one of '+', '.join(valid_plottypes)
    return

  # import the project communicator
  #from lasif.components.project import Project
  comm = Project(Lasif_path).comm

  # obtain information about station locations
  stations = comm.query.get_all_stations_for_event(event_name)
  nsta = len(stations)

  # obtain information about the windows
  wm_it1 = comm.windows.get(event_name, iteration1)
  wm_it2 = comm.windows.get(event_name, iteration2)
  
  # station list for windows
  window_channels_it1 = wm_it1.list()
  window_channels_it2 = wm_it2.list()
  
  if len(window_channels_it1) == 0:
    print "Event "+event_name+" contained no windows in iteration "+iteration+" in project "+Lasif_path
    return
  if len(window_channels_it2) == 0:
    print "Event "+event_name+" contained no windows in iteration "+iteration2+" in project "+Lasif_path
    return

  #- obtain station info and misfit info
  #  (I copied/adapted this from plot_stations_for_event)
  lngs = []
  lats = []
  lngs_nan = []
  lats_nan = []
  station_ids = []
  misfits = []
  visualisation_metric = []  #  can I immediately make ndarrays out of these? 
                             #+ Saves some converting later on when the faffing 
                             #+ with leaving list elements out is relevant (the isnan stuff)).
                             
                               # initialise the saving of extremal values
  extreme_pos = []
  extreme_neg = []

  if veryverbose:
    print 'Calculating misfits etc'
  #- for every station for that particular event, obtain misfits
  for stationID, value in stations.iteritems():  # Loop as dicts are unordered.
      
#        print stationID
      # get windows for station
      wins_it1 = wm_it1.get_windows_for_station(stationID)
      wins_it2 = wm_it2.get_windows_for_station(stationID)
      
      # get window info
      # (first determine whether there are windows at all (also fuck ugly))
      sz1 = np.size(wins_it1)  # sz will always be max 3: for each channel Z/N/E one
      sz2 = np.size(wins_it2)
      if sz1 > 0 and sz2 > 0:
      
        # go through both iterations to get the misfit change
        misf = {}; vizmet={}
        for it in [iteration1, iteration2]:
          
          if it == iteration1:
            sz = sz1; wins = wins_it1
          elif it == iteration2:
            sz = sz2; wins = wins_it2
          
          # get the misfits of all windows
          misfitsum = 0.; viz_metric_sum = 0.
          nwins = 0
          for ii in np.arange(sz):
            channelID = wins[ii].channel_id
            
            for jj in np.arange(len(wins[ii].windows)):
              nwins += 1
          
              # if plotting misfit...
              if plottype == 'misfits':
                misfit_val = wins[ii].windows[jj].misfit_value # non-normalised
                misfit_wt  = wins[ii].windows[jj].weight # window weight?
                misfit_weighted = misfit_val * misfit_wt 
                misfitsum += misfit_weighted
              
              
              # if plotting something else (phaseshift_integral, mean_delay, median_delay)
              else:
#                print wins[ii].windows[jj].misfit_details
                if plottype in wins[ii].windows[jj].misfit_details:
                  viz_metric = float(wins[ii].windows[jj].misfit_details[plottype])
                else:
                  print 'Warning: no '+plottype+' saved for channel '+channelID+', window {} !'.format(jj)
                  viz_metric = 0;
                viz_metric_sum += viz_metric      
            
                      
          # calculate _average_ of metric over all windows
          viz_metric_avg = viz_metric_sum/nwins
            
          misf[it] = misfitsum
          vizmet[it] = viz_metric_avg
            
        diff_misfit = misf[iteration2] - misf[iteration1]
        diff_viz_metric = abs(vizmet[iteration2]) - abs(vizmet[iteration1])
        
        # add misfit/phaseshift to array
        misfits.append(diff_misfit)
        visualisation_metric.append(diff_viz_metric)
        
        # add lonlats to arrays
        lngs.append(value["longitude"])
        lats.append(value["latitude"])
        station_ids.append(stationID)
        
        if not np.isnan(colour_max):
          if plottype=='misfits':
            if diff_misfit > colour_max:
              extreme_pos.append(stationID)
            elif diff_misfit < 0:
              extreme_neg.append(stationID)
          else:
            if diff_viz_metric > colour_max:
              extreme_pos.append(stationID)
            elif diff_viz_metric < -colour_max:
              extreme_neg.append(stationID)
        
      else:
        lngs_nan.append(value["longitude"])
        lats_nan.append(value["latitude"])
        #print "no windows for station "+stationID

  nwinsta = len(lngs)
  
  
  #- prepare the plotted shite:
  #   - put in a ndarray as otherwise indexing is bothersome
  #   - also make a color scale with 0 either as minimum or middle value
  import matplotlib.colors as colors
  from matplotlib import cm
  if plottype == "misfits":
    plotted_shite = np.asarray(misfits)
    collabel = "change in total misfit value (TF phase misfit, Fichtner et al. GJI 2008)"
    
  else:
    plotted_shite = np.asarray(visualisation_metric)
    collabel = "change in station averaged absolute {} (TF phase misfit, Fichtner et al. GJI 2008)".format(plottype.replace('_', ' '))
  
  if np.isnan(colour_max):
    maks = np.max(np.abs( plotted_shite[~np.isnan(plotted_shite)] ))
  else:
    maks = colour_max
  colnorm = colors.Normalize(vmin=-maks, vmax=maks)
  colmap = cm.RdBu_r
  colrange = [-maks, maks]
    
    
  print "colour maximum is "+str(maks)  
  
  #print plotted_shite
  
  
  
  ###################################################
  ## PLOTTING
  ###################################################
  
  iteration_name= iteration2 +' vs. '+ iteration1
  fig = plot_misfits_map(comm, iteration_name, event_name, plottype, 
                         lngs, lats, lngs_nan, lats_nan, 
                         plotted_shite, colmap, colnorm, colrange, collabel,
                         plot_windowless, plot_beachball, plot_histogram,
                         verbose=veryverbose)
  
  
  
  # add total misfit change to histogram
  if plot_histogram:
    (ax1, ax2, ax3) = fig.get_axes()
    total_value = plotted_shite.sum()
    iksleebel = ax2.get_xlabel()
    ax2.set_xlabel(iksleebel.replace('station averaged', 'change in station averaged absolute'))
#    ax2.set_xlabel('change in station averaged absolute '+plottype.replace('_', ' ')+' (values beyond '+str(colrange)+' binned in last bins!)')
    ax2.text(0.01, 0.85, 'Total change in '+plottype.replace('_', ' ')+' from iter '+iteration1+' to iter '+iteration2+': '+'{0:+.3f}'.format(total_value), transform=ax2.transAxes)

  # save & show
  #plt.tight_layout()
  if showdontsave:
    plt.ion()
    plt.show()
    plt.ioff()
  else:
    savename = 'misfit_change.'+iteration1+'-to-'+iteration2+'.'+event_name+'.'+plottype+'.avg'
    if plot_windowless:
      savename += '.with-winless'
    if plot_histogram:
      savename += '.with-hist'
    savename += '.png'
    savepath = os.path.join(savepath, savename)
    #print 'saving to: '+savepath
    plt.savefig(savepath)
    print 'saved figure to '+savepath
    plt.close()
    
  extreme_neg = sorted(extreme_neg)
  extreme_pos = sorted(extreme_pos)
  extreme_val = sorted(extreme_neg + [s + '  +' for s in extreme_pos])

  if verbose:
    print 'stations with extremely strong change in {} are:\n'.format(plottype.replace('_',' '))+'\n'.join(extreme_val)
  
  if output:
    return extreme_neg, extreme_pos
