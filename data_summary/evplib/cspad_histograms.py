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
            pylab.subplot(4,4,subplotnum)
            self.drawone(k)

    def drawone(self,k):
        pylab.bar(self.bins[:-1],self._histos[k],width=self.bins[1]-self.bins[0])
        pylab.xlim(self.bins[0],self.bins[-1])
        pylab.ylim(0.1,1.1*self._histos[k].max())



def splittile(arr):
    return arr[:,:194], arr[:,194:]

class cspadhisto(event_process.event_process):
    def __init__(self):
        self.bins = numpy.arange(-100,5000,200)
        self.histo = histo(64,self.bins)
        return

    def set_stuff(self,psana_src,psana_device,in_report=None,in_report_title=None):
        self.src         = psana.Source(psana_src)
        self.dev         = psana_device

    def event(self,evt):
        try:
            cspad = evt.get(self.dev, self.src)
        except Exception as e:
            return
        a = []
        for i in range(0,4):
            quad = cspad.quads(i)
            d    = quad.data()
            a.extend([ d[j] for j in range(0,8) ])

        for i, p in enumerate(a) :
            tile1, tile2 = splittile(p)
            self.histo.fill(2*i  ,tile1)
            self.histo.fill(2*i+1,tile2)

        return


if __name__ == "__main__":
    pylab.ion()
    ds = psana.DataSource('exp=CXI/cxig1515:run=173')
    psana.setConfigFile( 'cxicspad.cfg' )
    cs = cspadhisto()
    cs.set_stuff('DetInfo(CxiDs1.0:Cspad.0)' ,psana.CsPad.DataV2) # for this run use Ds2
    #cs.set_stuff('cspad_mod.CsPadCalib' ,psana.CsPad.DataV2) # for this run use Ds2
    events = ds.events()
    for i in xrange(100):
        evt = events.next()

    for i in xrange(30):
        cs.event(evt)

    #cs.histo.draw()
    cs.histo.drawone(0)

