#!/user/bin/env python3



from icecube.icetray import I3Tray,I3Units
from icecube import icetray, dataclasses, dataio
from icecube import tableio, hdfwriter



def getEvent(filePath,eventID,frameType=icetray.I3Frame.DAQ):
  "returns event with eventID from i3 file"
  for f in dataio.I3File(filePath):
    if f.Stop == frameType:
      if f["I3EventHeader"].event_id == eventID:
        break
  return f

def getGain(dom_status,omkey):
  return dom_status[omkey].dom_gain_type

def getGainGCD(GCD,omkey):
  f = dataio.I3File(GCD)
  dframe = f.pop_frame()
  while f.more() and not "I3DetectorStatus" in dframe:
      dframe = f.pop_frame()
  domStatus = dframe["I3DetectorStatus"].dom_status
  return getGain(domStatus,omkey)



def AddTotalTankHit(frame,pulseseriesList):
  '''calculates total SLC or HLC charges in tank pulses
  keys:[OfflineIceTopHLCTankPulses,OfflineIceTopSLCTankPulses,OfflineIceTopHLCVEMPulses,OfflineIceTopSLCVEMPulses]
  '''
  # print("Adding tank Hit")
  for pulseseries in pulseseriesList:
    NCh = 0
    psm = frame[pulseseries]
    if psm.__class__ == dataclasses.I3RecoPulseSeriesMapMask or psm.__class__ == dataclasses.I3RecoPulseSeriesMap:
      psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulseseries)
    for om,pulses in psm:
      for pulse in pulses:
        NCh +=1
        break
    channels = [omkey for omkey,ps in psm if len(ps)>0]
    # if NCh == 0:
    #   # print("psm",psm)
    frame[str(pulseseries)+"TotalHit"] = dataclasses.I3Double(NCh)