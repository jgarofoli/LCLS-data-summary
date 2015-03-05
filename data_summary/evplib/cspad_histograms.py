import os
import psana
import numpy
import logging
import event_process
import pylab
from mpi4py import MPI
from common import strtype

class result(object):
    """
    self describing, self reducing, self output generating
    generic result object
    """
    def __init__(self):
        return

class histo(object):
    def __init__(self,n,bins):
        self._histos = {}
        self._nentries = {}
        self._axes = {}
        self.bins = bins
        for i in xrange(n):
            self._histos[i] = numpy.histogram([],bins=self.bins)[0]
            self._nentries[i] = 0
        return

    def fill(self,i,arr):
        self._histos[i] += numpy.histogram(arr,bins=self.bins)[0]
        self._nentries[i] += 1
        return

    def draw(self):
        for k in sorted(self._histos):
            subplotnum=k%16
            if k % 16 == 0:
                pylab.figure()
            self._axes[k] = pylab.subplot(4,4,subplotnum+1)
            if subplotnum in [12,13,14,15]:
                pylab.xlabel('[{:0.0f}adc]'.format(self.bins[1]-self.bins[0]))
            if subplotnum in [0, 4, 8, 12]:
                pylab.ylabel('[N]')
            if subplotnum == 3:
                pylab.title('{:0.0f} frames'.format(self._nentries[k]))

            self.drawone(k)

    def drawone(self,k):
        pylab.plot(self.bins[:-1],self._histos[k])#,width=self.bins[1]-self.bins[0])
        pylab.xlim(self.bins[0],self.bins[-1])
        pylab.ylim(0.1,2*self._histos[k].max())
        pylab.gca().set_yscale('log')
        pylab.tick_params(axis='both', which='major', labelsize=6)
        pylab.tick_params(axis='both', which='minor', labelsize=6)
        pylab.gca().text(0.7,0.8,"tile {:0.0f}".format(k),fontsize=8,transform=pylab.gca().transAxes)
        #pylab.gca().set_xticklabels(   [int(i.get_text())/20. for i in pylab.gca().get_xticklabels() if len(i.get_text())>0] )



def splittile(arr):
    return arr[:,:194], arr[:,194:]

class cspadhisto(event_process.event_process):
    def __init__(self):
        self.bins = numpy.arange(-100,800,50)
        self.histo = histo(64,self.bins)
        return

    def set_stuff(self,psana_src,psana_device,in_report=None,in_report_title=None):
        self.src         = psana.Source(psana_src)
        self.dev         = psana_device

    def event(self,evt):
        try:
            self.cspad = evt.get(self.dev, self.src,'calibrated')
        except Exception as e:
            print e
            return 0

        if self.cspad is None:
            print "no cspad"
            return 0

        for i, p in enumerate(self.cspad) :
            tile1, tile2 = splittile(p)
            self.histo.fill(2*i  ,tile1)
            self.histo.fill(2*i+1,tile2)

        return 1

    def showcspad(self,ax=None):
        if not ax:
            ax = pylab.figure().gca()
        a = []
        for i in range(0,4):
            a.append( numpy.vstack( [cs.cspad[ i * j ] for j in range(0,8) ] ) )
        rawframe = numpy.hstack(a)
        ax.imshow(rawframe)

        


if __name__ == "__main__":
    psana.setConfigFile("cspad_ndarray.cfg") # must be set before defining DataSource
    pylab.ion()
    ds = psana.DataSource('exp=CXI/cxig1515:run=47') # 51, 52, 54
    #run=1
    #ds = psana.DataSource('exp=cxi86415:run=%d'%run)
    cls = ds.env().calibStore()
    cs = cspadhisto()
    cs.set_stuff('DetInfo(CxiDs1.0:Cspad.0)' ,psana.ndarray_float64_3) # for this run use Ds2
    #cs.set_stuff('cspad_mod.CsPadCalib' ,psana.CsPad.DataV2) # for this run use Ds2
    events = ds.events()

    cnt = 0
    #for i in xrange(1000):
    i = 0
    for evt in events:
        print i
        i += 1
        #evt = events.next()
        cnt += cs.event(evt)
        #if cnt >= 800:
        #    break

    cs.histo.draw()
    #cs.histo.drawone(0)

