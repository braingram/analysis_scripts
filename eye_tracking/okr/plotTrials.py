#!/usr/bin/env python

from pylab import *
from scipy.stats import linregress

from monkeyworks.data import MWKFile

timeSpan = (2091571591L, 2121555477L)
#timeSpan = (2967819517L, 2997803799L)

bestTrials = [107, 138]
goodTrials = [93, 94, 96, 106, 107, 134, 135, 138, 139, 143]

def read_okr_trial(timeSpan, dataFileName='/Users/graham/Data/20091028/okr_test_20091028_B.mwk', prePostTime=2000000):
    # get events for 1 sec prior to and 1 sec after stimulus
    timeSpan = (timeSpan[0] - prePostTime, timeSpan[1] + prePostTime)

    f = MWKFile(dataFileName)
    f.open()

    # data in okr_test_20091028_B.mwk suffer from the h/v swap bug
    hEvents = f.get_events(codes=['gaze_v',], time_range=timeSpan)
    vEvents = f.get_events(codes=['gaze_h',], time_range=timeSpan)
    
    aes = f.get_events(codes=['#announceStimulus',], time_range=timeSpan)
    stimProperties = {}
    for e in aes:
        if e.value['name'] == 'grating':
            stimProperties = e.value
    del aes
    
    f.close()
    
    def events_to_values(events):
        values = [e.value for e in events]
        times = [e.time for e in events]
        return times, values

    hTimes, hAngles = events_to_values(hEvents)
    vTimes, vAngles = events_to_values(vEvents)

    # resample
    times = zeros(len(hTimes) + len(vTimes), dtype=float64)
    angles = zeros((len(hTimes) + len(vTimes),2), dtype=float64)

    hIndex = 0
    vIndex = 0
    tIndex = 0
    wasH = False
    wasV = False
    while tIndex < len(times):
        if (hIndex >= len(hTimes)):
            newTime = vTimes[vIndex]
            newHAngle = hAngles[-1]
            newVAngle = vAngles[vIndex]
            vIndex += 1
        elif (vIndex >= len(vTimes)):
            newTime = hTimes[hIndex]
            newHAngle = hAngles[hIndex]
            newVAngle = vAngles[-1]
            hIndex += 1
        elif (hTimes[hIndex] < vTimes[vIndex]):
            if wasH == True: print "was H twice"
            wasH = True
            wasV = False
            # vtime is later than htime
            newTime = hTimes[hIndex]
            newHAngle = hAngles[hIndex]
            vM = (vAngles[vIndex] - vAngles[vIndex-1]) / (vTimes[vIndex] - vTimes[vIndex-1])
            newVAngle = vM * (newTime - vTimes[vIndex-1]) + vAngles[vIndex-1]
            #newVAngle = (vAngles[vIndex-1] + vAngles[vIndex]) / 2.0
            hIndex += 1
        else:
            if wasV == True: print "was V twice"
            wasV = True
            wasH = False
            # htime is later than vtime
            newTime = vTimes[vIndex]
            newVAngle = vAngles[vIndex]
            hM = (hAngles[hIndex] - hAngles[hIndex-1]) / (hTimes[hIndex] - hTimes[hIndex-1])
            newHAngle = hM * (newTime - hTimes[hIndex-1]) + hAngles[hIndex-1]
            #newHAngle = (hAngles[hIndex-1] + hAngles[hIndex]) / 2.0
            vIndex += 1
    
        times[tIndex] = newTime
        angles[tIndex] = newHAngle, newVAngle
        tIndex += 1

    # convert times from us to seconds
    times /= 1000000.0

    # remove blinks
    velocities = zeros(angles.shape, dtype=float64)
    # calculate velocities
    velocities[1:][:] = diff(angles, axis=0)
    dTimes = diff(times)
    velocities[1:,0] = velocities[1:,0] / dTimes
    velocities[1:,1] = velocities[1:,1] / dTimes
    velocity_threshold = 300.0
    vMask = abs(velocities) >= velocity_threshold
    errorSamples = vMask[:,0] + vMask[:,1]
    #print velocities.max()

    if any(errorSamples):
        errorIndices = errorSamples.nonzero()[0]
        print str(len(errorIndices))+" Errors found in file, attempting to fix..."
        i = 0
        while i < len(errorIndices):
            # find error window
            lastValid = errorIndices[i] - 1
            if i != len(errorIndices) -1:
                while (errorIndices[i+1] == errorIndices[i] + 1):
                    i += 1
                    if i == len(errorIndices) - 1:
                        break
            firstValid = errorIndices[i] + 1
            print lastValid, firstValid
            # if trying to look for a valid point out of bounds, loop around
            if firstValid >= len(times):
                firstValid = 0
            # interpolate points
            print "old points:", angles[lastValid+1:firstValid,:]
            hM = (angles[firstValid,0] - angles[lastValid,0])/(times[firstValid] - times[lastValid])
            vM = (angles[firstValid,1] - angles[lastValid,1])/(times[firstValid] - times[lastValid])
            dts = times[lastValid+1:firstValid] - times[lastValid]
            angles[lastValid+1:firstValid,0] = dts * hM + angles[lastValid,0]
            angles[lastValid+1:firstValid,1] = dts * vM + angles[lastValid,1]
            # angles[lastValid+1:firstValid,0] = linspace(angles[lastValid,0], angles[firstValid,0], \
            #             num=(firstValid-lastValid+1))[1:-1]
            # angles[lastValid+1:firstValid,1] = linspace(angles[lastValid,1], angles[firstValid,1],\
            #             num=(firstValid-lastValid+1))[1:-1]
            print "new points:", angles[lastValid+1:firstValid,:]
            i += 1
    
        # recalculate velocity
        velocities[1:][:] = diff(angles, axis=0)
        velocities[1:,0] = velocities[1:,0] / dTimes
        velocities[1:,1] = velocities[1:,1] / dTimes
    return times, angles, velocities, stimProperties

def get_trial_data(whiteList):
    ioff()
    f = open('trialTimes', 'r')
    trials = []
    for (i, l) in enumerate(f):
        if not i in whiteList:
            continue
        timeSpan = (int(l.split(' ')[0][1:-2]), int(l.split(' ')[1][:-3]))
        times, angles, velocities, stim = read_okr_trial(timeSpan)
        trials += [{'t': times, 'h': angles[:,0], 'v': angles[:,1],
                    'velH': velocities[:,0], 'velV': velocities[:,1],
                    'stim': stim, 'index': i},]
    return trials

def find_falls(values):
    i = 1
    falls = []
    start = None
    # find first fall
    while values[i] >= values[i-1]:
        i += 1
        if i == len(values):
            return []
    
    while i < len(values):
        if values[i] >= values[i-1]:
            if start != None:
                falls += [{'start': start, 'length': i - 1 - start,
                            'amp': values[i-1] - values[start]},]
                start = None
            i += 1
        else:
            if start == None:
                start = i - 1
            i += 1
    return falls

def find_rises(values):
    i = 1
    rises = []
    start = None
    # find first fall
    while values[i] <= values[i-1]:
        i += 1
        if i == len(values):
            return []
    
    while i < len(values):
        if values[i] <= values[i-1]:
            if start != None:
                rises += [{'start': start, 'length': i - 1 - start,
                            'amp': values[i-1] - values[start]},]
                start = None
            i += 1
        else:
            if start == None:
                start = i - 1
            i += 1
    return rises

def measure_okr(h, t, falls, minPPoints=30, minSPoints=3, minP=0.1, figNum=None):
    pursuitVel = []
    saccadeVel = []
    i = 1
    if figNum != None:
        figure(figNum)
    while i < len(falls):
        pStart = falls[i-1]['start']+falls[i-1]['length']+2
        pEnd = falls[i]['start']
        sStart = falls[i]['start']
        sEnd = falls[i]['start']+falls[i]['length']
        if pEnd - pStart < minPPoints or sEnd - sStart < minSPoints:
            i += 1
            continue
        pr = linregress(t[pStart:pEnd],h[pStart:pEnd])
        sr = linregress(t[sStart:sEnd],h[sStart:sEnd])
        if (pr[3] <= minP and sr[3] <= minP):
            saccadeVel += [sr[0],]
            pursuitVel += [pr[0],]
            if figNum != None:
                ts = array([t[sStart],t[sEnd]])
                ys = ts * sr[0] + sr[1]
                plot(ts,ys,c='g',linewidth=2)
                ts = array([t[pStart],t[pEnd]])
                ys = ts * pr[0] + pr[1]
                plot(ts,ys,c='r',linewidth=2)
        i += 1
    if figNum != None:
        plot(t,h,c='k')
    return pursuitVel, saccadeVel

def analyze_okr(trial,figNum=None):
    if trial['stim']['direction'] == 0:
        # look for rises
        moves = find_rises(trial['h'])
    else:
        # look for falls
        moves = find_falls(trial['h'])
    moves = [m for m in moves if (abs(m['amp']) > 1.0)]
    if len(moves) == 0:
        return "No rises/falls found"
    pV, sV = measure_okr(trial['h'], trial['t'], moves, figNum=figNum)
    speed = trial['stim']['speed']
    closedLoopGain = abs(median(pV))/speed
    openLoopGain = abs(median(sV))/speed
    yl = ylim()
    vlines([trial['t'][m['start']] for m in moves], yl[0], yl[1])
    ylim(yl)
    print "%i\t%.2f\t%.2f" % (trial['index'], closedLoopGain, openLoopGain)
    return closedLoopGain, openLoopGain, pV, sV

if __name__ == '__main__':
    ioff()
    f = open('trialTimes','r')
    imgDir = 'images/'
    i = 0
    #whiteList = bestTrials
    whiteList = goodTrials
    #whiteList = []
    for l in f:
        # add possible whiteList
        if len(whiteList) != 0:
            if not i in whiteList:
                i += 1
                continue
        # if i < 111:
        #     i += 1
        #     continue
        timeSpan = (int(l.split(' ')[0][1:-2]), int(l.split(' ')[1][:-3]))
        #print 'reading trial...'
        times, angles, velocities, stim = read_okr_trial(timeSpan)
        #print stim
        # def chop(N,alpha=0.05):
        #             st = arange(0, len(t) - N, N)
        #             slopes = []
        #             for i in xrange(len(t)-(N+1)):
        #                 r = lr(t[i:i+N+1],a[i:i+N+1,0])
        #                 if (r[3] < alpha): slopes.append(r[0])
        #             return slopes
        #         
        # plot okr trial
        figure(1)
        suptitle(str(timeSpan)+' D:'+str(stim['direction'])+' SF:'+str(stim['frequency'])+' Speed:'+str(stim['speed']))
        subplot(221); plot(times - times[0], angles[:,0]); title('Horizontal Gaze Angle'); ylim((-40.,0.))
        xlabel('Time (seconds)'); ylabel('Gaze Angle (degrees)')
        vlines((2.0,(times[-1]-times[0]-2.0)),ylim()[0],ylim()[1])
        subplot(222); plot(times - times[0], velocities[:,0]); title('Horizontal Eye Velocity'); ylim((-250.,250.))
        xlabel('Time (seconds)'); ylabel('Eye velocity (degrees/second)')
        vlines((2.0,(times[-1]-times[0]-2.0)),ylim()[0],ylim()[1])
        subplot(223); plot(times - times[0], angles[:,1]); title('Vertical Gaze Angle'); ylim((-40.,0.))
        xlabel('Time (seconds)'); ylabel('Gaze Angle (degrees)')
        vlines((2.0,(times[-1]-times[0]-2.0)),ylim()[0],ylim()[1])
        subplot(224); plot(times - times[0], velocities[:,1]); title('Vertical Eye Velocity'); ylim((-250.,250.))
        xlabel('Time (seconds)'); ylabel('Eye velocity (degrees/second)')
        vlines((2.0,(times[-1]-times[0]-2.0)),ylim()[0],ylim()[1])
        # aesthetic adjustments
        subplots_adjust(wspace=0.34,hspace=0.41)
        # save figure
        savefig(imgDir+str(i)+'.jpg', dpi=300)
        close(1)
        i += 1
        #break
    f.close()
    
