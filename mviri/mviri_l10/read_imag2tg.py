"""
this file contains functions reading the imag2tg file
"""

import struct
from collections import namedtuple
import numpy as np
import glob
import mmap
import os
import sys
try:# if available: load matplotlib
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt
  import matplotlib.mlab as mlab
  import matplotlib
except ImportError:
  print "NO PLOTTING WILL BE POSSIBLE"

global recordsize

recordsize=17296 #records are described further down in structures()


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]


def l15_to_L10_name(rect2lpname,s_datestring):
  """
  This function searches the corresponding marf2tg filename 
  for the given rect2lp filename. Note that the marf2tg file 
  then first has to be converted to imag2tg
  Input:
   -rect2lpname, rect2lp filename including full path
  Output:
   -imag2tgname, imag2tg filename including full path
  """
  prepend="/"#"/tcenas/proj/eraclim/data/imag2tg" #necessary while converted imag2tg are not yet available on NAS
  folders=rect2lpname.split("/")
  outfolder=""
  for folder in folders[:-1]:
    if "level1" in folder:
      outfolder="/".join((outfolder,"level0"))
    elif "HR-MFG15" in folder:
      outfolder="/".join((outfolder,"IMAG2TG"))
    else:
      #outfolder="/".join((outfolder,folder))
      outfolder=os.path.join(outfolder,folder)
    #print outfolder
  idx=folders[-1].index(s_datestring)
  outfolder="/".join((outfolder,folders[-1][idx+6:idx+8])) #METEOSAT7-MVIRI-MTP15-NA-NA-20070731233000
  outfolder="/".join((outfolder,"*"+folders[-1][idx:idx+12]+"*"))
  if not prepend=="":
    outfolder=prepend+outfolder
  outfolder.replace("//","/")
  imag2tgname=glob.glob(outfolder)
  #print imag2tgname[0]
  try:
    return imag2tgname[0]
  except IndexError:
    print "ERROR: no fitting L1.0 file found using wildcard: "
    print outfolder
    exit()


def decode_telem(lininfo,head,outfile="out.csv",mode="psi"):
  """
  This function extracts and decodes the telemetry information from the 
  line headers.
  Input:
  -lininfo, line headers/trailers as read by function image()
  -head, image headers as read by function header()
  -outfile, filename to store the decoded TM frames; obsolete if not "s" in mode
  -mode, string to determine what to do with the TM info:
    containing "p" for plot, "s" for storing, i for read from image, 
    m for read from mmap and/or "r" for return to caller. 
    either i or m have to be specified!!
  Output:
  -outtab,only returned if mode contains "r"
  Example1:
      >>>import readerL1 as rd
      >>>f=open(filename, "rb")
      >>>f,head = rd.header(f)
      >>>f,lininfo,pix,pixir,pixwv = rd.image(f,head)
      >>>f.close()
      >>>buffer=rd.decode_telem (lininfo,head,mode="ri")
  Example2:
      >>>with open(iname, "rb") as f:
      >>>  f,head0   = r0.header(f)
      >>>  lininfo0  = r0.telem_mmap2(f)
      >>>  telem=r0.decode_telem(lininfo0,head0,mode="rm")
  """
  
#READ TELEMETRY OFFSETS
  CONFIGPP=os.environ['CONFIGPP']
  teletab=np.genfromtxt(CONFIGPP+"/Telemetry_offsets.csv",delimiter=",",dtype="str",skip_header=1)
  telem_descr=dict(tuple(teletab[:,0:2]))
  recN=2500#3029-20
  #convert teletab rows into outtab columns
  outtab=np.zeros((recN/32+10),dtype={'names':np.append(['index',"jday","HH","MM","SS","QA","missing"],teletab[:,0]),'formats':np.repeat('i4',len(teletab[:,0])+7)})
  teleframe= np.empty((0), int)
#GO THROUGH ALL LINES
  ld=-1
  idx=0
  missing=0
  tmfr_old=0
  qual=0
  for li in range(recN):
    #print li
    if "i" in mode:
      teleHead = np.array(getattr(lininfo[li], 'byTelemHead')).astype("ubyte")
    elif "m" in mode:
      teleHead = np.array(lininfo[li,0:20]).astype("ubyte")
#CHECK IF TM BELONGS TO SLOT  
    TTAG=decode_timetag(teleHead)
    # Calculate the slot in day from the timetag
    TTAG_SLOTINDAY = 1 + (2 * TTAG[1]) + (TTAG[2]/30)
    HEAD_SLOTINDAY = getattr(head[0], 'iSlot_no')[0]
    if TTAG_SLOTINDAY == HEAD_SLOTINDAY:
      if ld==-1:
        ld=li #store first valid line
      l=li-ld
#IF BELONGS TO SLOT: TEST AND PROCESS
    if ld>=0:
      J_SAT = teleHead[3-1]
      J_DAY = getattr(head[0], 'iJulian_day')[0]
      #print teleHead
      if "i" in mode:
        telbuff=getattr(lininfo[li], 'byTelemetry')
      elif "m" in mode:
        telbuff=lininfo[li,20:52]
      tmfr=telbuff[2]#framenumber
      #save frame in bits
      telerec = np.unpackbits(np.array(telbuff).astype("ubyte"))
      #read tm frame quality
      if not qual==0:
        qual=int("".join(telerec[0:16].astype(str)),2)
      #eventually fill missing frames
      if li>0:
        tz=0
        while (tmfr>tmfr_old+1+tz) and not (tmfr_old+1+tz>32):
          missing=missing+1
          qual=-1
          teleframe=np.append(teleframe,np.unpackbits(np.repeat(0,32).astype("ubyte")))
          tz=tz+1
      #check quality
      frameQA=framecheck(J_SAT,tmfr,TTAG_SLOTINDAY,HEAD_SLOTINDAY,TTAG[0],J_DAY)
      if frameQA==0:#checks tm frame contents
#READ AND STORE
        #if current is a new tm frame:
        if (li>0) and ((tmfr==0) or ((tmfr==1) and (tmfr_old >= 29)) or ((tmfr==2) and (tmfr_old >= 30)) or((tmfr==3) and (tmfr_old == 31))):
          if len(teleframe)!=10496:
            qual=-3
  #CONTENT DEFINITION FROM TELETAB
          #do deciphering
          #split outtab according to teletab rows and fill
          for line in teletab:
            st=(int(line[2])*32+int(line[3]))*8+int(line[4])
            #print line[0]+" ("+line[1]+") starts at line:"+line[2]+" offset: "+str(st)
            en=st+int(line[5])
            #print np.shape(teleframe)
            bitstring=teleframe[st:en]
            #print line[0]+" ("+line[1]+") is: "+str(bitstring)
            intstring=np.packbits(bitstring)
            #print line[0]+" ("+line[1]+") is: "+str(intstring)
            #print len(teleframe)
  #STORE
            #test length and store actual variable
            try:
              outtab[line[0]][idx]=intstring[0]
            except IndexError:
              outtab[line[0]][idx]=-1
              if qual<0:#because TM frame not complete: modify error code
                qual=qual-2
              else:
                qual=-2 
          #store some metadata
          outtab['index'][idx]=l
          outtab['jday'][idx]=J_DAY
          outtab['HH'][idx]=TTAG[1]
          outtab['MM'][idx]=TTAG[2]
          outtab['SS'][idx]=TTAG[3]
          outtab['missing'][idx]=missing
          outtab["QA"][idx]=qual
#RESTART NEXT TM FRAME
          #start new
          idx=idx+1
          missing=0
          teleframe=telerec
          qual=0
#FILL CURRENT FRAME
        else:
          teleframe=np.append(teleframe,telerec)
        tmfr_old=tmfr
          #print np.shape(teleframe)
          
      else:
        print "ERROR: wrong TM frame or bad quality!"
        print " occured at line "+str(l)+"("+str(l%32)+") with tm frame: "+str(tmfr) + " and qual: "+str(qual)
  if "p" in mode:
    y="y"
    while y=="y":
      vari=raw_input("whattcha wanna plot?")
      plt.plot(outtab["index"][outtab["QA"]==0],outtab[vari][outtab["QA"]==0],'.')
      plt.show()
      y=raw_input("Whanna plot more? (y/n)")
  if "s" in mode:
    np.savetxt(outfile,outtab,header=";".join(outtab.dtype.names),delimiter=";",comments='')
  if "r" in mode:
    return outtab,telem_descr

def framecheck(J_SAT,J_FRAMENO,TTAG_SLOTINDAY,J_SLOTINDAY,TTAG_DAY,J_DAY):
  """
  Helper function for decode_telem()
  """
  if J_SAT == 0 or J_SAT == 255 :
    q=1
  elif J_FRAMENO < 0 or J_FRAMENO > 31 :
    q=2
  elif ( TTAG_SLOTINDAY != J_SLOTINDAY ) and ( TTAG_SLOTINDAY != J_SLOTINDAY - 1 and not (TTAG_SLOTINDAY == 48 and J_SLOTINDAY == 1) ):
    q=3
  elif ( TTAG_DAY != J_DAY ) and ( TTAG_DAY != J_DAY - 1  and not ( TTAG_DAY >= 365 and J_DAY == 1) ):
    q=4
  else:
    q=0
  return q

def decode_timetag(HKDATA):
  """
  Helper function for decode_telem()
  """
  #! Decode the timetag - not all parts actually needed here
  ITEMP1 = HKDATA[13-1]/16               # x100 DAY
  ITEMP2 = HKDATA[13-1] - (ITEMP1 * 16)  # x10 day
  ITEMP3 = HKDATA[14-1]/16               # x1 day
  TTAG_DAY = (ITEMP1 * 100) + (ITEMP2 * 10) + ITEMP3

  ITEMP1 = HKDATA[14-1] - (ITEMP3 * 16)  # x10 hour
  ITEMP2 = HKDATA[15-1]/16               # x1 hour
  TTAG_HOUR  = (ITEMP1 * 10) + ITEMP2

  ITEMP1 = HKDATA[15-1] - (ITEMP2 * 16)  # x10 min
  ITEMP2 = HKDATA[16-1]/16               # x1 min
  TTAG_MIN = (ITEMP1 * 10) + ITEMP2

  ITEMP1 = HKDATA[16-1] - (ITEMP2 * 16)  # x10 sec
  ITEMP2 = HKDATA[17-1]/16               # x1 sec
  TTAG_SEC = (ITEMP1 * 10) + ITEMP2

  ITEMP1 = HKDATA[17-1] - (ITEMP2 * 16)  # x100 msec
  ITEMP2 = HKDATA[18-1]/16               # x10 msec
  ITEMP3 = HKDATA[18-1] - (ITEMP2 * 16)  # x1 msec
  TTAG_MSEC = (ITEMP1 * 100) + (ITEMP2 * 10) + ITEMP3
  
  return TTAG_DAY,TTAG_HOUR,TTAG_MIN,TTAG_SEC,TTAG_MSEC

def telem_mmap2(f):#for quickly stripping the tm frames
  #faster
  headers=39
  image_offset=recordsize*headers
  line_offset=17208#[12,5012]
  fl=0
  ll=2500
  tm_frame=[]
  for l in range(fl,ll):
    offset=image_offset+l*recordsize+line_offset
    data = np.memmap(f, dtype='byte', mode='r', offset=offset,shape=(52))
    tm_frame.append(data.copy())
  return np.array(tm_frame)


def space_corner_mmap(f,debug=0):#for quickly stripping the spacecorners
  #slower than mmap2
  headers=39
  image_offset=recordsize*headers
  sensors=[0,1]
  line_offsets=[160,5160]#[12,5012]
  lines =np.hstack((range(0,200),range(2300,2500)))
  pixels=np.hstack((range(0,400),range(4600,5000)))
  space=[]
  for sensor in sensors:
    spacevalues=[]
    if debug==1:
      offsetz=[]
    for l in lines:
      for p in pixels:
	offset=image_offset+l*recordsize+line_offsets[sensor]+p
	mm = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
	data = mm[offset:offset+1]
	spacevalues.append(struct.unpack('b', data)[0])
	if debug==1:
	  offsetz.append(offset)
    space.append(spacevalues)
  if debug==1:
    return np.array(space),np.array(offsetz)
  else:
    return np.array(space)

def space_corner_mmap2(f,debug=0):#for quickly stripping the spacecorners
  #faster
  headers=39
  image_offset=recordsize*headers
  sensors=[0,1]
  line_offsets=[160,5160]#[12,5012]
  lb_l=np.array([0,0,2300,2300])
  ub_l=lb_l+200
  lb_p=np.array([0,4600,0,4600])
  ub_p=lb_p+400
  space=[]
  for sensor in sensors:
    spacevalues=[]
    if debug==1:
      offsetz=[]
    for corner in np.arange(0,4):
      for l in range(lb_l[corner],ub_l[corner]):
	p=lb_p[corner]
	offset=image_offset+l*recordsize+line_offsets[sensor]+p
	data = np.memmap(f, dtype='byte', mode='r', offset=offset,shape=(400))
	spacevalues=np.hstack((spacevalues,data))
	if debug==1:
	  offsets=np.add(np.arange(lb_p[corner],ub_p[corner]),image_offset+l*recordsize+line_offsets[sensor])
	  offsetz=np.hstack((offsetz,offsets))

    space.append(spacevalues)
  if debug==1:
    return np.array(space),np.array(offsetz)
  else:
    return np.array(space)

def space_corner_mmap3(f,debug=0):#for quickly stripping the spacecorners
  #should be fastest
  headers=39
  image_offset=recordsize*headers
  sensors=[0,1]
  line_offsets=[160,5160]#[12,5012]
  lines =np.hstack((range(0,200),range(2100,2500)))
  pixels=np.hstack((range(0,400),range(4600,5000)))
  p,l=np.meshgrid(pixels,lines)
  space=[]
  for sensor in sensors:
    offsets=image_offset+l*recordsize+line_offsets[sensor]+p
    data = np.memmap(f, dtype='byte', mode='r', offset=np.ndarray.flatten(offsets),shape=np.shape(np.ndarray.flatten(offsets)))
    space.append(data)
  return np.array(space)
  
def header(f,start=1,end=39):
  #open file
    records=range(start,end+1)#records=[1,2,3]
    
    recs=[]
    for r in records:

       #handle records>39
      rq=r
      if r>39 and r<=3069:
	rq=40
      
      #get information on record structure    
      names,dims,rec_fmt=structures(rq)
      recID	='r'+str(r)
      struc 	= namedtuple(recID, names, verbose=False)
      
      #read and unpack record
      rec_bin = f.read(recordsize)
      #print str(rec_fmt[0]),r
      record=struct.unpack(str(rec_fmt[0]),rec_bin)
      
      #go through fields
      g=0
      lis=[]
      for i in range(0,len(names)):
	
	currname=names[i]
	currtype=currname[:1]
	buf=[]
	for j in range(dims[i]):
	  #buf.append(vax(record[g+j],currtype))
	  buf.append(record[g+j])
	  
	#collect content in one list
	lis.append(buf)
	g=g+dims[i]

      
      #assign content to named tuple
      struc=struc(*lis)
      
      #store record
      recs.append(struc)
    
    return f,recs  

def image(f,head,start=40,end=3069):

    vistargets=['shVisSPixel','shVisNPixel']
    irtargets =['shIRPixel']
    wvtargets =['shWVPixel']
    pixeldata=['shVisSPixel','shVisNPixel','shIRPixel','shWVPixel']
  
    #open file
    records=range(start,end+1)#records=[1,2,3]
    
    recs=[]
    visim=[]
    irim=[]
    wvim=[]
    for r in records:
	
       #handle records>39
      rq=r
      if r>39 and r<=3069:
	rq=40
      names,dims,rec_fmt=structures(rq)
      recID	='r'+str(r)
      struc 	= namedtuple(recID, names[:4]+names[8:], verbose=False)
      
      #read and unpack record
      rec_bin = f.read(recordsize)
      record=struct.unpack(str(rec_fmt[0]),rec_bin)
      
      #debug
      #print "\nrecord: "+str(r)+"\n"
      #print struct.calcsize(rec_fmt[0])#print size of the structure
      #print len(record)
      #print len(names)
      early=['M2','P2']
      sat=getattr(head[0], 'szSat_code')[0]

      #go through fields
      g=0
      lis=[]
      for i in range(0,len(names)): 
	currname=names[i]
	currtype=currname[:1]
	
	#debug
	#print "name: "+str(names[i])
	#print "dims: "+str(dims[i])
	#print "orig[0]: "+str(record[(g)])
	#print "conv[0]: "+str(vax(record[(g)],currtype))
	#print currtype
	
	#fill buffer to capture dimensional fields
	buf=[]
	for j in range(dims[i]):
	  #debug
	  #print j
	  #print "orig: "+str(record[g+j])
	  #print "conv: "+str(vax(record[g+j],currtype))
	  #print sat , early
	  if sat in early:
	    try:
	      buf.append(int(format(record[g+j], 'b').zfill(8)[:6],2))#
	    except ValueError:
	      #print "Value error: "+str(record[g+j])
	      buf.append(np.nan)
	  else:
	    try:
	      buf.append(int(record[g+j]))
	    except ValueError:
	      #print "Value error: "+str(record[g+j])
	      buf.append(np.nan)
	#store line headers and trailers in list...
	if not currname in pixeldata:
	  #print i,currname,pixeldata
	  lis.append(buf)
	#...and store pixel values in image
	elif currname in pixeldata:
	  if currname in vistargets:
	    visim.append(buf)
	  if currname in irtargets:
	    irim.append(buf)
	  if currname in wvtargets:
	    wvim.append(buf)
	g=g+dims[i]

      
      #assign content to named tuple
      struc=struc(*lis)
      #debug
      #print ', '.join(['{0}={1}'.format(k, getattr(struc, k)) for k in struc._fields])

      #store record
      recs.append(struc)
    print getattr(recs[2],"iJulianSlot")
    return f,recs,visim,irim,wvim


def structures(z):
  import re
  #define names
  str_r=[['szOriginalFormat','iYear','iJulian_day','iSlot_no','szSpare2','szSat_code','szSpare3'],
         ['szSpare'],
         ['szSpare'],
         ['iJulian_slot_no','shJulian_day','shSlot_no','shNominal_image_time',\
          'shSlot_type','shImage_quality_flag','szSat_id',\
          'shSub_sat_point_displacement_in_pixels','shSub_sat_point_displacement_in_lines',\
          'szVis_missing_line_table','szVisn_missing_line_table',\
          'szIr_missing_line_table','szWv_missing_line_table',\
          'shNo_missing_lines_replaced','shNo_black_lines_13','shNo_black_lines_23','shNo_black_lines_33',\
          'szWefax_annotation28','szImage_quality_annotation12',\
          'shSpectral_contents3','shDeformation_used',\
          'szBb_ir_count_for_space_view6','szStd_Bb_ir_count_for_space_view3',\
          'szBb_wv_count_for_space_view6','szStd_Bb_wv_count_for_space_view3',\
          'szBb_cold_temp5','szBb_warm_temp5','szBb_ir_count_nominal6','szStd_Bb_ir_count_nominal3',\
          'szBb_wv_count_nominal6','szStd_Bb_wv_count_nominal3','szBb_calibration_timestamp5',\
          'szMpef_abs_ir_cal5','szMpef_ir_space_count3','szMpef_ir_calibration_timestamp5',\
          'szMpef_abs_wv_cal5','szMpef_wv_space_count3','szMpef_wv_calibration_timestamp5',\
          'szSpares115','szAll_channel_gains8','szSpares24',\
          'dsRight_ascension_att_south','dsDeclination_att_south','dsRight_ascension_att_north','dsDeclination_att_north',\
          'dsRefined_attitude_xyz','dsRight_ascension_declination_of_mean_att2',\
          'iNo_slots_with_refined_attitude','flSpin_dur_minus_nominal_spin_dur',\
          'szEclipse_operation','szDecontamination','szManoeuvre','szView','szIr1_on','szIr2_on','szWv1_on','szWv2_on','szVis1_on','szVis2_on','szVis3_on','szVis4_on',\
          'szMpef_spares20','szSpares416','szImage_status16','shLine_pixel_orientation12',\
          'dsSat_earcth_centre_distance','dsOrbit_offset_at_southern_horizon',\
          'dsOrbit_offset_at_northern_horizon3','flMax_deformation_diff_x_inside_column',\
          'flMax_deformation_diff_y_inside_line','flMax_deformation_diff_x_inside_line',\
          'flMax_deformation_diff_y_inside_column','szImage_conditions4',\
          'shMin_count_in_histogram','shMax_count_in_histogram4','flMean_vis_1_4','flMean_vis_2_3','flSnr_in_space_corners',\
          'iNo_lines_in_snr_calc','flSnr_eastern_part','flSnr_western_part','flMean_noise_count_eastern_part','flMean_noise_count_western_part',\
          'shMax_space_count_eastern_part','shMax_space_count_western_part','szSpares5','szMessage',\
          'shNominal_sat_longitude','szSpares','flNominal_sat_longitude_in_degrees',\
          'iCode_for_sub_sat_point','iNo_landmarks','iNo_vis_cloudfree_landmarks','iNo_ir_cloudfree_landmarks',\
          'iGQASource','iHeaderVersionNumber','flNo_vis_landmarks_with_corr','flNo_ir_landmarks_with_corr',\
          'flAbs_std_visible_landmarks','flAbs_max_visible_landmarks','iNo_landmarks_with_max_deviation',\
          'flAbs_std_ir_landmarks','flAbs_max_ir_landmarks',\
          'flAbs_landmark_results','flNo_vis_landmarks_with_curr',\
          'flNo_ir_landmarks_with_curr','flRel_std_visible_landmarks',\
          'flRel_max_visible_landmarks','iNo_landmarks_with_rel_max_deviation',\
          'flRel_std_ir_landmarks','flRel_max_ir_landmarks','flRel_landmark_results',\
          'flVis_equalisation_coef','flDif_direct_in_pixel','flDif_direct_in_line',\
          'iRectification_flag','shGQA_check_flag','iImage_processing_status',\
          'szAlignPadding','flSpaceCornersVISS','flSpaceCornersVISN',\
          'flSpaceCornersIR','flSpaceCornersWV','flStdSpaceCornersVISS','flStdSpaceCornersVISN',\
          'flStdSpaceCornersIR','flStdSpaceCornersWV','szSpares9'],
          ['szSpare'],#5
          ['szSpare'],#6
          ['szSpare1','dsNominal_sat_longitude','szSpare2','dsOrbitCoordinatesFixedEarthImageStart','dsOrbitCoordinatesFixedEarthImageEnd''szSparez'],#7 
          ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], #15
          ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], ['szSpare'], #25
          ['szSpare'], ['szSpare'], #27
          ['iJulianSlot','flRTDeformMatr_X_1_35'], #28
          ['iJulianSlot','flRTDeformMatr_X_36_70'], #29
          ['iJulianSlot','flRTDeformMatr_X_71_110'], #30
          ['iJulianSlot','flRTDeformMatr_Y_1_35'], #31
          ['iJulianSlot','flRTDeformMatr_Y_36_70'], #32
          ['iJulianSlot','flRTDeformMatr_Y_71_110'], #33
          ['iJulianSlot','flBaDeformMatr_X_1_35'], #34
          ['iJulianSlot','flBaDeformMatr_X_36_70'], #35
          ['iJulianSlot','flBaDeformMatr_X_71_110'], #36
          ['iJulianSlot','flBaDeformMatr_Y_1_35'], #37
          ['iJulianSlot','flBaDeformMatr_Y_36_70'], #38
          ['iJulianSlot','flBaDeformMatr_Y_71_110'], #39
       ['iJulianSlot','shLineNr','shSILineNr','szSpares1',\
          'shVisSPixel','shVisNPixel','shIRPixel','shWVPixel',\
          'szSpares2','byTelemHead','byTelemetry','szSparessz']]

  dim_r= [[8,1,1,1,52,1,17222],
          [17296],
          [17296],#4152
          [1,1,1,1,1,1,1,1,1,316,316,316,316,1,3,3,3,1,12,3,1,6,3,6,3,\
            5,5,6,3,6,3,5,5,3,5,5,3,5,15,8,4,1,1,1,1,3,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,\
            20,16,16,12,1,3,3,1,1,1,1,4,4,4,1,1,4,1,4,4,4,4,4,4,88,1,1,2,1,1,1,\
            1,1,1,1,1,1,2,2,2,2,2,1024,1,1,2,2,2,2,2,1024,1,4,4,1,1,1,2,4,\
            4,4,4,4,4,4,4,6042],
          [17296],#5
          [17296],#6
          [11232,1,880,6,6,5080],#7 
          [17296], [17296], [17296], [17296], [17296], [17296], [17296], [17296], #15
          [17296], [17296], [17296], [17296], [17296], [17296], [17296], [17296], [17296], [17296], #25
          [17296], [17296], #27
          [1,105*35], #28
          [1,105*35], #29
          [1,105*35], #30
          [1,105*35], #31
          [1,105*35], #32
          [1,105*35], #33
          [1,105*35], #34
          [1,105*35], #35
          [1,105*35], #36
          [1,105*35], #37
          [1,105*35], #38
          [1,105*35], #39
          [1,1,1,152,5000,5000,2500,2500,2048,20,32,36]]#len:17296

  rec_fmt= [['<8ciii52c2s17222c'],
            ['<17296c'],
            ['<17296c'],
            ['<ihhhhh2shh316c316c316c316ch3h3h3h28s12b3hh6c3c6c3c\
              5c5c6c3c6c3c5c5c3c5c5c3c5c15c8c4cdddd3d2difbbbbbbbb\
              bbbb20c16c16b12hd3d3dffff4b4h4hff4fi4f4f4f4f4h4h88c\
              800sh2cf4siiiiifff2f2i2f2f1024fff2f2f2i2f2f1024ff4f\
              4fihi2c4f4f4f4f4f4f4f4f3d26c6042c'],
            ['<17296c'],#5
            ['<17296c'],#6
            ['<11232cd880c6d6d5080c'],#7 
            ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], #15
            ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], ['<17296c'], #25
            ['<17296c'], ['<17296c'], #27
            ['<i4323f'], #28
            ['<i4323f'], #29
            ['<i4323f'], #30
            ['<i4323f'], #31
            ['<i4323f'], #32
            ['<i4323f'], #33
            ['<i4323f'], #34
            ['<i4323f'], #35
            ['<i4323f'], #36
            ['<i4323f'], #37
            ['<i4323f'], #38
            ['<i4323f'], #39
            ['ihh152c5000B5000B2500B2500B2048c20B32B36c']]#40
  
  

  #replace non-alphanumeric values
  str_o=str_r[z-1]
  i=0
  for stri in str_o:
    str_o[i]=re.sub(r'\W+', '',str_o[i])
    i=i+1
  #return
  return (str_o,dim_r[z-1],rec_fmt[z-1])