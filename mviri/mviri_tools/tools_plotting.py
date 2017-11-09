"""
this file contains functions for plotting
"""

import matplotlib.pyplot as plt
import numpy as np

def printimage(pix,tit,maxi=255,save=0):
    """
    Simple plotting of an MVIRI image
    Inputs:
    pix   array with pixel values
    tit   title/name of the image
    """
    pix=np.array(pix)
    pix[pix>maxi]=maxi
    im=plt.imshow(np.uint8(np.fliplr(pix)),origin='lower', interpolation='none')
    if save == 0:
      plt.show()
    else:
      print "   to: ../tmp/"+tit+".png"
      plt.savefig("../tmp/"+tit+".png")
      plt.close()