import os
import psana
import numpy
import logging
import event_process
import pylab
from mpi4py import MPI
from common import strtype
import toolbox

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
        self._figs = []
        self._auxfigs = []
        self.bins = bins
        self._aux = {}  
        self._aux_index = []  
        for i in xrange(n):
            self._histos[i]   = numpy.histogram([],bins=self.bins)[0]
            self._nentries[i] = numpy.array([0])

            self._aux[i] = { 'index': [], 'sum': [], 'mean': [], 'std': [], 'median': [] }

        return

    def fill(self,i,arr):
        # b) sum, mean, std and median for each subpanel
        # 
        thishist = numpy.histogram(arr,bins=self.bins) 
        self._aux[i]['sum'].append( arr.sum() )
        self._aux[i]['mean'].append( arr.mean() )
        self._aux[i]['std'].append( arr.std() )
        self._aux[i]['median'].append( numpy.median(arr) )
        self._aux_index.append( None )
        self._histos[i] += thishist[0]
        self._nentries[i][0] += 1
        return

    def draw(self):
        for k in sorted(self._histos):
            subplotnum=k%16
            if k % 16 == 0:
                self._figs.append(pylab.figure())
            self._axes[k] = pylab.subplot(4,4,subplotnum+1)
            if subplotnum in [12,13,14,15]:
                pylab.xlabel('[{:0.0f}adc]'.format(self.bins[1]-self.bins[0]))
            if subplotnum in [0, 4, 8, 12]:
                pylab.ylabel('[N]')
            if subplotnum == 3:
                pylab.title('{:0.0f} frames'.format(self._nentries[k][0]))

            self.drawone(k)

    def drawone(self,k):
        x,y = toolbox.mksteps(self.bins, self._histos[k]+0.1,startval=0.1) 
        pylab.fill_between( x, 0.1, y[:-1] )
        #pylab.plot(self.bins[:-1],self._histos[k])#,width=self.bins[1]-self.bins[0])
        pylab.xlim(self.bins[0],self.bins[-1])
        pylab.ylim(0.1,2*self._histos[k].max())
        pylab.gca().set_yscale('log')
        pylab.tick_params(axis='both', which='major', labelsize=6)
        pylab.tick_params(axis='both', which='minor', labelsize=6)
        pylab.gca().text(0.65,0.8,"tile {:0.0f}".format(k),fontsize=8,transform=pylab.gca().transAxes)
        #pylab.gca().set_xticklabels(   [int(i.get_text())/20. for i in pylab.gca().get_xticklabels() if len(i.get_text())>0] )

    def draw2dhisto(self):
        self._auxfigs.append(pylab.figure())
        self.axis = self._auxfigs[-1].add_subplot(111)
        axis = self.axis
        pylab.imshow(  numpy.log10(numpy.vstack( numpy.array( [self._histos[i] for i in sorted(self._histos)] ) )+1.) , aspect='auto', interpolation='none' )
        pylab.xlabel('[arb]')
        pylab.ylabel('subpanel')
        pylab.title('2dhistogram')
        pylab.colorbar()
        #ax = axis.get_xticks(), axis.get_xticklabels()
        #ra = ax[0][-1]-ax[0][0]
        #br = self.bins[-1] - self.bins[0] + (self.bins[1]-self.bins[0])
        #self.newlabels = []
        #for i in xrange(len(ax[0])):
        #    if len(ax[1][i].get_text()) != 0:
        #        self.newlabels.append( "{:0.0f}".format((float(i)/len(ax[0])*br)+self.bins[0]) )
        #    else :
        #        self.newlabels.append('')
        #axis.set_xticklabels(self.newlabels)

    def drawaux(self,var):
        self._auxfigs.append(pylab.figure())
        pylab.imshow( numpy.transpose( numpy.vstack( numpy.array( [self._aux[i][var] for i in sorted(self._aux)] ) ) ) , aspect='auto' , interpolation='none')
        pylab.title(var)
        pylab.xlabel('subpanel')
        pylab.ylabel('event')
        pylab.colorbar()
        return

    def drawallaux(self):
        for v in sorted(self._aux[0]):
            if v == 'index':
                continue
            self.drawaux(v)




def splittile(arr):
    return arr[:,:194], arr[:,194:]

def product(*it):
    tot = 1
    for i in it:
        tot *= i
    return tot

class cspadhisto(event_process.event_process):
    def __init__(self,**kwargs):
        subpanels = kwargs.get('subpanels',64)
        self.bins      = kwargs.get('bins', numpy.arange(-100,800,5))
        self.logger                      = logging.getLogger(__name__+'.cspadhist')
        self.output = event_process.event_process_output()
        self.output['in_report']         = None
        self.output['in_report_title']   = None
        self.pixels = product( *(32,185,388) )
        if self.pixels % subpanels == 0: # should be gotten from first cspad
            self.subpanels = subpanels # must break up the cspad into and even number of subpanels
            self.panelsize = self.pixels/self.subpanels
        else:
            self.logger.error('subpanels doesn\'t subdivide pixels evenly') # more info
            self.subpanels = 64

        self.histo = histo(self.subpanels,self.bins) # get away from the hardcoded 64 here
        return

    def set_stuff(self,psana_src,psana_device,in_report=None,in_report_title=None):
        self.src         = psana.Source(psana_src)
        self.dev         = psana_device
        self.output['in_report']         = in_report
        self.output['in_report_title']   = in_report_title

    def event(self,evt):
        try:
                              # type (not device)
            self.cspad = evt.get(self.dev, self.src,'calibrated')
        except Exception as e:
            print e
            return False

        if self.cspad is None:
            self.logger.info( "no cspad" )
            return False

        #self.flat = self.cspad.flatten() # .reshape((1,product(*self.cspad.shape))) # what is the arrangement of these subpanels??
        self.flat = numpy.array( [ self.cspad[i,:,:].flatten() for i in range(self.cspad.shape[0])] ).flatten() # would like this step to be FASTER
        for i in xrange(self.subpanels) :
            self.histo.fill(i  ,self.flat[i*self.panelsize:(i+1)*self.panelsize])

        return True

    def showcspad(self,ax=None):
        if not ax:
            ax = pylab.figure().gca()
        a = []
        for i in range(0,4):
            a.append( numpy.vstack( [cs.cspad[ i * j ] for j in range(0,8) ] ) )
        rawframe = numpy.hstack(a)
        ax.imshow(rawframe, interpolation='none')
        pylab.colorbar()

    def endJob(self):
        self.logger.info('mpi reducing cspad histograms')
        self.merged_histograms = histo(self.subpanels,self.bins)

        for k in sorted(self.merged_histograms._histos):
            self.parent.comm.Reduce(self.histo._histos[k],   self.merged_histograms._histos[k],   op=MPI.SUM, root=self.reducer_rank)
            self.parent.comm.Reduce(self.histo._nentries[k], self.merged_histograms._nentries[k], op=MPI.SUM, root=self.reducer_rank)
            # reduce the aux info
        for i in sorted(self.histo._aux):
            for k in sorted(self.histo._aux[i]):
                data = self.parent.comm.gather(self.histo._aux[i][k], root=self.reducer_rank)
                for d in data:
                    self.merged_histograms._aux[i][k].extend(d)
        data = self.parent.comm.gather(self.histo._aux_index,root=self.reducer_rank)
        for d in data:
            self.merged_histograms._aux_index.extend(d)


        if self.parent.rank == self.reducer_rank:
            self.logger.info("cspad histograms reducing")
            self.output['text'].append("some text")

            self.merged_histograms.draw()
            self.merged_histograms.drawallaux()
            self.merged_histograms.draw2dhisto()

            self.output['figures'] = {}
            for i, fig in enumerate(self.merged_histograms._figs) :
                self.output['figures']['histos_{:0.0f}'.format(i)] = {}
                fig.savefig(os.path.join( self.parent.output_dir, 'figure_cspad-histo-{:0.0f}.png'.format(i) ), dpi=600)
                self.output['figures']['histos_{:0.0f}'.format(i)]['png'] = os.path.join( self.parent.output_dir, 'figure_cspad-histo-{:0.0f}.png'.format(i) )

            for i, fig in enumerate(self.merged_histograms._auxfigs) :
                title = fig.get_axes()[0].get_title()
                self.output['figures']['aux_{:}'.format(title)] = {}
                fig.savefig(os.path.join( self.parent.output_dir, 'figure_cspad-aux-{:}.png'.format(title) ), dpi=600)
                self.output['figures']['aux_{:}'.format(title)]['png'] = os.path.join( self.parent.output_dir, 'figure_cspad-aux-{:}.png'.format(title) )

            self.parent.output.append(self.output)

        return

        


if __name__ == "__main__":
    psana.setConfigFile("cspad_ndarray.cfg") # must be set before defining DataSource
    pylab.ion()

    if True :
        ds = psana.DataSource('exp=CXI/cxig1515:run=47') # 51, 52, 54
        cs = cspadhisto(subpanels=64,bins=numpy.arange(-100,400,5))
        cs.set_stuff('DetInfo(CxiDs1.0:Cspad.0)' ,psana.ndarray_float64_3) # for this run use Ds2
    else :
        run=6
        ds = psana.DataSource('exp=cxi86415:run=%d'%run)
        cs = cspadhisto(subpanels=64,bins=numpy.arange(-100,400,5))
        cs.set_stuff('DetInfo(CxiDs2.0:Cspad.0)' ,psana.ndarray_float64_3) # for this run use Ds2
    #cls = ds.env().calibStore()
    #cs.set_stuff('cspad_mod.CsPadCalib' ,psana.CsPad.DataV2) # for this run use Ds2
    events = ds.events()

    cnt = 0
    i = 0
    Ntot = 100
    for evt in events:
        print i
        i += 1
        #evt = events.next()
        cnt += cs.event(evt)
        if cnt >= Ntot and True :
            break

    #cs.histo.draw()
    cs.histo.drawallaux()
    cs.histo.draw2dhisto()
    #cs.histo.drawone(0)

