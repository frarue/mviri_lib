import mviri.mviri_l10 as mviri_l10
print "readers imported"
import mviri.mviri_tools as mviri_tools
print "tools imported"

print "\nExample: \n/DSNNAS/Repro_Temp/mviri/level0/IMAG2TG/data/MET2/1982/08/08/METEOSAT2-MVIRI-MTP10-NA-NA-19820808113000.000000000Z"
iname=raw_input("\nPlease provide sample filepath+filename as in example above: \n")


#read imag2tg++++++++++
print "reading ..."
with open(iname, "rb") as f:
  f,head0                           =   mviri_l10.header(f)
  f, lininfo0, pixvis, pixir, pixwv =   mviri_l10.image(f,head0)   #reads line headers/trailers AND pixel values
  #lininfo0                         =   mviri_l10.telem_mmap2(f) #reads only line headers/trailers -->faster!!

#plot VIS image
print "plotting ..."
#mviri_tools.printimage(pixvis,"VIS - "+ iname)

#extract telemetry+++++
print "extracting telemetry ..."
telem,telem_descr=mviri_l10.decode_telem(lininfo0,head0,mode="rm")

#read space corners++++
print "extracting space corners ..."
with open(iname, "rb") as f:
  space = mviri_l10.space_corner_mmap2(f)                         #quickly strips space corners



