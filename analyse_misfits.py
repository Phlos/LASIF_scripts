        # scripts to analyse misfit

import matplotlib.pyplot as plt
import numpy as np
import os
from lasif.components.project import Project
from lasif import LASIFNotFoundError
from lasif import LASIFAdjointSourceCalculationError

def calc_misfit_for_channel(comm, event_name, iteration_name, channel_name):

  '''
  calculate the misfit for a specific channel
  '''
  
  # obtain information about the windows
  wm = comm.windows.get(event_name, iteration_name)
  
  it = comm.iterations.get(iteration_name)
  event_weight = it.events[event_name]["event_weight"]
  
  wins = wm.get(channel_name)
    
  station_weight = it.events[event_name]["stations"][".".join(channel_name.split(".")[:2])]["station_weight"]
  
  channel_misfit = 0
  total_channel_weight = 0
  
  for win in wins.windows:
  
    try:
      misfit = win.misfit_value
    except LASIFAdjointSourceCalculationError:
      continue
    except LASIFNotFoundError as e:
      print str(e)
      continue
    
    channel_misfit += misfit * win.weight
    total_channel_weight += win.weight
  
  # rare, but sometimes all wins for a certain channel fail the calculation
  if total_channel_weight == 0:
    print 'total channel weight is zero'
    
  # make sure the misfits are consistent with adj src calculation
  channel_misfit *= event_weight * station_weight / total_channel_weight
  print "channel {1} misfit is {0:.5f}".format(channel_misfit, channel_name)
  
  
  

def calc_misfit_per_event(comm, iteration_name, event_name):
  
  # loop over windows to calculate total misfit
  
  # obtain information about the windows
  wm = comm.windows.get(event_name, iteration_name)
  # obtain information about all stations for this event
  #stations = comm.query.get_all_stations_for_event(event_name)
  
  it = comm.iterations.get(iteration_name)
  event_weight = it.events[event_name]["event_weight"]
  
  
  channels = wm.list()
  if len(channels) == 0:
    print 'Warning! event '+event_name+' has no windows for iteration '+iteration_name+' ! '
    return np.nan
  
  event_misfit = 0
  for channel in channels:
    
    wins = wm.get(channel)
    
    station_weight = it.events[event_name]["stations"][".".join(channel.split(".")[:2])]["station_weight"]
    
    channel_misfit = 0
    total_channel_weight = 0
    
    for win in wins.windows:
    
      try:
        misfit = win.misfit_value
      except LASIFAdjointSourceCalculationError:
        continue
      except LASIFNotFoundError as e:
        print str(e)
        continue
      
      channel_misfit += misfit * win.weight
      total_channel_weight += win.weight
    
    # rare, but sometimes all wins for a certain channel fail the calculation
    if total_channel_weight == 0:
      continue
      
    # make sure the misfits are consistent with adj src calculation
    channel_misfit *= event_weight * station_weight / total_channel_weight
    #print "channel {1} misfit is {0:.2f}".format(channel_misfit, channel)
    event_misfit += channel_misfit
    
  return event_misfit
  
  '''
  # station list for windows
  window_stations = wm.list()

  if len(window_stations) == 0:
    print "Event "+event_name+" contained no windows in iteration "+iteration+" in project "+Lasif_path
    return 0.
    
  #- for every station for that particular event, obtain misfits
  misfit_sum = 0.
  for stationID, value in stations.iteritems():  # Loop as dicts are unordered.
      
    # get windows for station
    wins = wm.get_windows_for_station(stationID)
    
    # get window info
    

    # get the misfits of all windows
    # determine number of windows
    sz = np.size(wins)
    misfit_sta = 0.
    for ii in np.arange(sz):
    
      misfit_val = wins[ii].windows[0].misfit_value # non-normalised
      misfit_wt  = wins[ii].windows[0].weight # window weight?
      misfit_weighted = misfit_val * misfit_wt 
      misfit_sta += misfit_weighted
      
    misfit_sum += misfit_sta
  '''  
  return misfit_sum

def calc_misfit(Lasif_path='./EMed_full/', iteration_name='0_new', save=False, verbose=True):
  
  '''
  For a given Lasif project and iteration name, calculate the total misfit
  '''
  
  # import the project communicator
  comm = Project(Lasif_path).comm
  
  it = comm.iterations.get(iteration_name)
  
  # for each event, calculate event misfit (loop over windows like in plot_misfit_map etc)
  total_misfit = 0.
  dict_event_misfits = {}
  for event_name in sorted(it.events.keys()):
    if verbose:
      print ' ...event '+event_name
    event_misfit = calc_misfit_per_event(comm, iteration_name, event_name)
    dict_event_misfits[event_name] = event_misfit
    total_misfit += event_misfit
    
  if save:
    outfile = 'iter_misfit.'+iteration_name+'.npy'
    print 'saving file to '+outfile
    np.save(outfile, (total_misfit, dict_event_misfits) )
    
  return total_misfit, dict_event_misfits

def calc_misfit_development_dict(misfit_dict, event_list, iter_list):
  
  '''
  Dictionary transform tool. Super dumb.
  
  Make a dict of misfits with order misfit_dict[event_name][iteration_name]
  from one with order misfit_dict[iteration_name][event_name]
  
  input:
  misfit_dict:  stored as misfit_dict[iter_name][event_name]
  '''
  

  dict_misfits_per_event = {}
  for event_name in event_list:
  
    misf_per_event = {}
    for iteration_name in iter_list:
      # the index [1] is because index [0] stores the total iteration misfit
      if event_name in misfit_dict[iteration_name][1].keys():
        misf_per_event[iteration_name] = misfit_dict[iteration_name][1][event_name]
    
    if any(misf_per_event):
      # add to dictionary if there are any iterations with this event
      dict_misfits_per_event[event_name] = misf_per_event
          
  return dict_misfits_per_event 

#def add_sisyphus_misfit_log(ax1, filename, normalise):
  '''
  This function does not work well
  '''

  #misfit_log = [float(x.split(' ')[-1].strip()) for x in open(filename).readlines()]
  
  #iters=[]
  #n_full_iters = len(misfit_log)/4
  #for ii in range(n_full_iters):
    #for jj in range(5)[1:]:
      #iters.append(ii+0.1*jj)
  
  ##print 'iters len {}, misfit log len {}'.format(len(iters), len(misfit_log))
  
  ## plotting
  #print 'value of misfit_log[0] is {}'.format(misfit_log[0])
  #print 'value of misfit_log[1] is {}'.format(misfit_log[1])
  #if normalise:
    #misfit_log = [x / misfit_log[0] for x in misfit_log]
  #ax1.plot(iters, misfit_log, label='SISYPHUS misfits incl. test steps')
  #ax1.plot(iters[3::4], misfit_log[3::4], marker='o', label='SISYPHUS misfits of _4 models')
  
  #return iters, misfit_log

def plot_misfit_development(Lasif_path='EMed_full', iter_list=['10_new', '19_4'], 
                            misfit_dict_in={}, highlight_event=[],
                            normalise=True, sort_alpha_by='init_misfit', verbose=True,
                            save=False, savepath='./plots/'):
  
  '''
  Plot the misfit developement as a function of iteration
  for the supplied list of iterations (or 'ALL')
  
  INPUT:
  lasif_path      the path to the Lasif project
  iter_list = []  a list that contains the iterations to be 
                  compared. They will be shown in the order supplied
  misfit_dict={}  a list of pre-calculated misfits per iter. If it
                  is empty, the misfits will be calculated
  normalise       should all misfits be normalised by their initial
                  value.
  sort_alpha_by   how should the event misfit lines be shaded?
  verbose         How much output is given
  save            True or False; whether the figure should be plotted
                  or saved.
  savepath:       location where plot is saved ('./plots/' by default)
  
  OUTPUT:
  dict_misfits_per_event :    a dictionary of misfits per event
                              ordered misfit{event}{iter}
  misfit_dict :               same, but ordered misfit{iter}{event}
  {a figure} :                a figure with misfit development per 
                              event as a function of iteration.
  
  '''
  
  ##############################
  ### some preparation
  ##############################
  
  # make sure the function doesn't accidentally work on external vars.
  misfit_dict = misfit_dict_in.copy()
  
  # import the project communicator & get event list
  comm = Project(Lasif_path).comm
  event_list = sorted(comm.events.list())
  
  
  ############################################
  ### Loading/calculating event misfits
  ############################################
  
  # loop over iters to check if misfit exists
  for it in iter_list:
    misfit_npy_file_name='iter_misfit.'+it+'.npy'
    if it in misfit_dict.keys():
      print 'Misfit for iteration '+it+' already exists'
      continue
    elif os.path.exists(misfit_npy_file_name):
      print 'Loading misfit for iteration '+it+' from file'
      misfit_dict[it] = np.load(misfit_npy_file_name)
    else:
      print 'Calculating misfit for iteration '+it
      misfit_dict[it] = calc_misfit(Lasif_path, it, save=True, verbose=verbose)
  
  
  ############################################
  ### Make new dict with misf[event][iter] 
  ###         instead of misf[iter][event]
  ############################################
  
  
  #  (This method is objectively rather stupid... 
  #+ I wish I had a table or something)
  
  # invert the dictionary
  #if not any(dict_misfits_per_event):
  dict_misfits_per_event = calc_misfit_development_dict(misfit_dict, event_list, iter_list)
  
  # loop over iterations to obtain iteration misfit list
  iter_misfit_total=[]
  for iter_name in iter_list:
    
    if iter_name in misfit_dict.keys():
      iter_misfit_total.append(misfit_dict[iter_name][0])
    else:
      iter_misfit_total.append(np.nan)
  
  # loop over events to gather the development of misfits
  #dict_misfits_per_event = {}
  # also save some other event properties for plotting sorting
  misf1 = {}; maxmisf1 = 0; 
  event_mags={}; maxmag = 0; minmag=15
  evmag_normd = {}; maxmag_normd = 0
  no_misfits = []
  
  for event_name in event_list:
    
    if event_name in dict_misfits_per_event.keys():
      misf_per_event = dict_misfits_per_event[event_name]
    else:
      print 'event '+event_name+' has no misfits for any of these iterations.'
      no_misfits.append(event_name)
      
  # make a list of events that are in at least one of the iterations  
  event_list_reduced = [x for x in event_list if x not in no_misfits]
  
  # calculate some properties
  for event_name in event_list_reduced:
    # first misfit determination
    for iter_name in iter_list:
      if event_name in misfit_dict[iter_name][1].keys():
        misf1[event_name] = misfit_dict[iter_name][1][event_name]
        break
    maxmisf1 = max(misf1[event_name], maxmisf1)
    # event magnitude
    event_mags[event_name] = comm.events.get(event_name)['magnitude']
    maxmag = max(event_mags[event_name], maxmag)
    minmag = min(event_mags[event_name], minmag)
    # number of stations
    nsta = len(comm.query.get_all_stations_for_event(event_name))
    # normalised event magnitude
    evmag_normd[event_name] = event_mags[event_name] / nsta
    maxmag_normd = max(maxmag_normd, evmag_normd[event_name])
    
  #------------------------
  # PLOTTING
  #------------------------
  
  if save:
    plt.ioff()
  
  
  fig, ax1 = plt.subplots(1, figsize=(12.,10.))
  
  # plot in figure
  plotline = {}; 
  firstline=True; 
  for event_name in event_list_reduced:
      
    # transform dict to lists (for plotting)
    y=[]
    for itname in iter_list:
      if itname in dict_misfits_per_event[event_name]:
        y.append(dict_misfits_per_event[event_name][itname])
      else:
        print 'Warning! event '+event_name+' has no misfit for iter '+itname
        y.append(np.nan)
    x = np.arange(len(iter_list)) + 0.5
    
    # plot real or normalised misfit
    if normalise:
      plotline[event_name] = ax1.plot(x,y/misf1[event_name], color='black')
      labely=('normalised misfit')
    else:
      plotline[event_name] = ax1.plot(x,y, color='black')
      labely=('misfit')
      
    # only give one event a label (not 83 black lines)
    if firstline:
      plotline[event_name][0].set_label('misfits per event, shaded by '+sort_alpha_by)
      firstline=False
      
    # shade the lines according to sort_alpha_by
    if sort_alpha_by=='init_misfit':
      alpha_value = misf1[event_name]/maxmisf1
    elif sort_alpha_by=='event_mag':
      alpha_value = (event_mags[event_name]-minmag)/(maxmag-minmag)
      alpha_value = max(alpha_value, 0.1)
      #print 'alpha for '+event_name+' is {0:.2f} because mag here is {1:.1f}'.format(alpha_value, event_mags[event_name])
    elif sort_alpha_by=='event_mag_by_nsta':
      alpha_value = (evmag_normd[event_name] / maxmag_normd )
    plotline[event_name][0].set_alpha(alpha_value)
    
    # give the highlighted event special looks
    if event_name in highlight_event:
      plotline[event_name][0].set_color('green')
      plotline[event_name][0].set_linewidth(2)
      plotline[event_name][0].set_label('event '+event_name)
    
  
  # Plot total misfit development
  if normalise:
    ax1.plot(x,iter_misfit_total/iter_misfit_total[0], color='red', linewidth=2, label='total misfit')
    eilim = ax1.get_ylim()
    ax1.set_ylim(max(eilim[0], 0.7), min(eilim[1], max(iter_misfit_total/iter_misfit_total[0]) * 1.1))
  else:
    ax1.plot(x,iter_misfit_total, color='red', linewidth=2, label='total misfit')
    eilim = ax1.get_ylim()
    ax1.set_ylim(max(eilim[0], 0.7*min(iter_misfit_total)), min(eilim[1], max(iter_misfit_total) * 1.1))
    
  
  # Prettify the plot
  ax1.set_xticks(x)
  ax1.set_xticklabels(iter_list)
  ax1.set_xlim(0,len(iter_list))
  ax1.set_xlabel('iteration name')
  ax1.set_ylabel(labely)
  ax1.set_title('Misfit development across iterations')
  plt.legend(loc='upper left')
  
  
  
  
  # actual drawing
  if save:
    if normalise:
      normorno='.normalised'
    else:
      normorno=''
    if len(highlight_event)>0:
      hili='.highlight-'+str(highlight_event)
    else:
      hili=''
    savename='misfit_development.'+iter_list[0]+'-to-'+iter_list[-1]+normorno+'.alpha_sorted_'+sort_alpha_by+hili+'.png'
    savepath=os.path.join(savepath,savename)
    print 'saving figure to '+savepath
    plt.savefig(savepath)
    plt.close()
  else:
    plt.ion()
    plt.show()
    plt.ioff()
    
  
  return dict_misfits_per_event, misfit_dict

def plot_misfit_log_sisyphus(filename, iter_list, 
                             save=False, savepath='./plots/', inversion_name='inv_20170731'):
  
  '''
  Plot the misfit development based on Alexey's misfit
  '''
      
  misfit_log = [x.split(' ')[-1].strip() for x in open(filename).readlines()]
  
  iters=[]
  n_full_iters = len(misfit_log)/4
  for ii in range(n_full_iters):
    for jj in range(5)[1:]:
      iters.append(ii+0.1*jj)
  
  print 'iters len {}, misfit log len {}'.format(len(iters), len(misfit_log))
  # plotting
  
  fig, ax1 = plt.subplots(1, figsize=(12.,10.))
  
  ax1.plot(iters, misfit_log, label='misfits incl. test steps')
  ax1.plot(iters[3::4], misfit_log[3::4], marker='o', label='misfits of _4 models')
  
  # Prettify the plot
  ax1.set_xlabel('iteration')
  ax1.set_ylabel('total misfit')
  ax1.set_title('Misfit development from SISYPHUS misfit.log for '+inversion_name)
  ax1.set_xticks(iters[3::4])
  ax1.set_xticklabels(iter_list)
  plt.legend()
  
  # Actual displaying
  if save:
    savename=inversion_name+'.misfit_development.sisyphus.png'
    savepath=os.path.join(savepath,savename)
    print 'saving figure to '+savepath
    plt.savefig(savepath)
    plt.close()
  else:
    plt.ion()
    plt.show()
    plt.ioff()
    
def get_good_bad_events(misfit_dict, iter1, iter2, goodorbad=None):
  
  '''
  Print out the events that are good or bad (misfit increase, or misfit decrease)
  '''
  
  good_ones = []; bad_ones=[]
  for event_name in misfit_dict[iter1][1].keys():
    m_init = misfit_dict[iter1][1][event_name]
    m_final= misfit_dict[iter2][1][event_name]
    if m_final < m_init:
      good_ones.append(event_name)
    elif m_final > m_init:
      bad_ones.append(event_name)

  if goodorbad=='bad':
    for event_name in bad_ones:
      m_init = misfit_dict[iter1][1][event_name]
      m_final= misfit_dict[iter2][1][event_name]
      incr = m_final - m_init
      rel_incr = (m_final - m_init) / m_init
      print 'event {0:<9} is bad: '.format(event_name)+ \
            'misfit went {0:6.2f} --> {1:6.2f}'.format(m_init, m_final) +\
            '({0:+5.1f}  = {1:+5.1f}%)'.format(incr, 100.*rel_incr)
  elif goodorbad=='good':
    for event_name in good_ones:
      m_init = misfit_dict[iter1][1][event_name]
      m_final= misfit_dict[iter2][1][event_name]
      incr = m_final - m_init
      rel_incr = (m_final - m_init) / m_init
      print 'event {0:<9} is good: '.format(event_name)+ \
            'misfit went {0:6.2f} --> {1:6.2f}'.format(m_init, m_final) +\
            '({0:+5.1f}  = {1:5.1f}%)'.format(incr, 100.*rel_incr)
          
  return good_ones, bad_ones
